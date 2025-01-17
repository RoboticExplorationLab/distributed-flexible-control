function [timeVec, XVec, uAll, zSol, fval] = ...
    solveMultiArcManeuverMS(problemStruct, options)
% SOLVEMULTIARCMANEUVERMS
%
% Multiple-shooting solver for a multi-arc rest-to-rest (or general) maneuver.
%
% All inputs are packed into 'problemStruct', which should have fields:
%   .xInit   (nx x 1) initial state
%   .xTarget (nx x 1) target/final state
%   .dynFun(t, X, u)          : user-supplied continuous dynamics => dX/dt
%   .discreteDynFun(X, u, dt) : user-supplied discrete dynamics => X_next
%   .arcConstraintFun(tArray, XArray) => [c, ceq]
%   .costFun(tAll, XAll) => scalar
%   .uList       (1 x N) array of piecewise-constant controls
%   .guessT      (N x 1) initial guess for arc durations
%   .guessX      (nx x N) guess states at each arc boundary (X2..X_{N+1})
%
% The user can either:
%   - Provide .dynFun for continuous integration (old default), OR
%   - Provide .discreteDynFun for discrete updates.
%   If both are given, .discreteDynFun takes precedence.
%
% 'options' is the fmincon options struct or [] for defaults.
%
% Outputs:
%   timeVec : row vector of all time samples across arcs
%   XVec    : nx x length(timeVec) array of states
%   uAll    : 1 x length(timeVec) array of piecewise-constant controls
%   zSol    : final decision variable vector
%   fval    : final cost

    % If options is empty, set defaults
    if isempty(options)
        options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
        options.MaxFunctionEvaluations = 1e6;
        options.ScaleProblem = 'obj-and-constr';
        options.ConstraintTolerance = 1e-2 * options.ConstraintTolerance;
        options.MaxIterations = 4000;
    end

    % Extract data from the struct
    xInit            = problemStruct.xInit;
    xTarget          = problemStruct.xTarget;
    arcConstraintFun = problemStruct.arcConstraintFun;
    costFun          = problemStruct.costFun;
    uList            = problemStruct.uList;
    guessT           = problemStruct.guessT;
    guessX           = problemStruct.guessX;

    %                We'll pick continuous vs discrete in subfunctions
    if isfield(problemStruct,'dynFun')
        contDynFun = problemStruct.dynFun;  % continuous
    else
        contDynFun = []; % not provided
    end
    
    if isfield(problemStruct,'discreteDynFun') && ~isempty(problemStruct.discreteDynFun)
        discDynFun = problemStruct.discreteDynFun; % discrete
    else
        discDynFun = []; % not provided
    end

    if isfield(problemStruct,'finalConstraintFun') && ~isempty(problemStruct.finalConstraintFun)
        finalConstraintFun = problemStruct.finalConstraintFun;
    else
        % Default final constraint = Xf - xTarget = 0
        finalConstraintFun = @(Xf) deal([], Xf - xTarget(:));
    end

    N  = length(uList);      % number of arcs
    nx = length(xInit);      % dimension of the state

    % Pack initial guess for fmincon
    z0 = packZ(guessT, guessX);

    % Bounds for time >= 0, no upper bound; states unbounded
    lbT = zeros(N,1);
    ubT = problemStruct.maxT*ones(N,1);
    lbX = -inf(nx*N,1);
    ubX = +inf(nx*N,1);
    lb  = [lbT; lbX];
    ub  = [ubT; ubX];
    
    % Example linear constraint that sum of all times <= 300
    A = [ones(1,N), zeros(1,nx*N)];
    b = problemStruct.maxT;

    % Create the optimization problem
    if problemStruct.multistart
        problem = createOptimProblem('fmincon', ...
            'objective',  @(z) objWrapper(z, xInit, uList, costFun, contDynFun, discDynFun), ...
            'x0',         z0, ...
            'Aineq',      A, ...
            'bineq',      b, ...
            'lb',         lb, ...
            'ub',         ub, ...
            'nonlcon',    @(z) constrWrapper(z, xInit, uList, ...
                                            arcConstraintFun, finalConstraintFun, ...
                                            contDynFun, discDynFun), ...
            'options',    options);
    
        % MultiStart
        ms = MultiStart('Display','iter','UseParallel',true,'MaxTime',3600);
        numStartPoints = problemStruct.NMStart;
        [zSol, fval] = run(ms, problem, numStartPoints);
    else
        [zSol, fval] = fmincon( ...
                     @(z) objWrapper(z, xInit, uList, costFun, contDynFun, discDynFun), ...
                     z0, ...
                     A, b, ...     % no linear constraints
                     [], [], ...
                     lb, ub, ...
                     @(z) constrWrapper(z, xInit, uList, ...
                                            arcConstraintFun, finalConstraintFun, ...
                                            contDynFun, discDynFun), ...
                     options);
    end

    % Reconstruct the full trajectory from the solution
    [timeVec, XVec] = buildFullTrajectory(zSol, xInit, uList, contDynFun, discDynFun);

    % Build piecewise-constant control profile
    uAll = buildControlProfile(timeVec, zSol, uList, nx);
end


%% =====================================================================
%% Objective function wrapper
function Jval = objWrapper(z, xInit, uList, costFun, contDynFun, discDynFun)
    [tAll, XAll] = buildFullTrajectory(z, xInit, uList, contDynFun, discDynFun);
    Jval         = costFun(tAll, XAll);
end

%% =====================================================================
%% Nonlinear constraints wrapper
function [c, ceq] = constrWrapper(z, xInit, uList, ...
                                  arcConstraintFun, finalConstraintFun, ...
                                  contDynFun, discDynFun)
    [tVec, Xbndry] = unpackZ(z, length(uList), length(xInit));
    N = length(uList);

    cList   = [];
    ceqList = [];

    Xi = xInit;  % boundary state for arc i (X1)
    for i = 1:N
        dt = tVec(i);
        ui = uList(i);

        % Integrate from Xi over [0, dt] with control ui
        [tArc, XArc] = simulateArc(Xi, dt, ui, contDynFun, discDynFun);

        % End-state from integration
        Xend = XArc(:, end);

        % Next boundary from decision vector
        XiPlus1 = Xbndry(:, i);

        % Enforce continuity: Xi+1 - Xend = 0
        ceqList = [ceqList; (XiPlus1 - Xend)];

        % Arc constraints from user
        [cArc, ceqArc] = arcConstraintFun(tArc, XArc);
        cList   = [cList;   cArc(:)];
        ceqList = [ceqList; ceqArc(:)];

        % Move to next arc boundary
        Xi = XiPlus1;
    end

    % Final-state constraint: X_{N+1} = xTarget
    % ceqFinal = Xi - xTarget(:);
    [cEnd, ceqEnd] = finalConstraintFun(Xi);
    cList   = [cList;  cEnd(:)];
    ceqList = [ceqList; ceqEnd];

    c   = cList;
    ceq = ceqList;
end

%% =====================================================================
%% Single-arc integrator
function [tSol, XSol] = simulateArc(x0, dt, u, contDynFun, discDynFun)
    % If dt <= 0, trivial arc
    if dt <= 0
        tSol = 0;
        XSol = x0;
        return;
    end

    if ~isempty(discDynFun)
        %%% <-- DISCRETE DYNAMICS
        % Example: one-step update for the entire "arc"
        % tSol = [0, dt], XSol = [x0, x_next]
        Xnext = discDynFun(x0, u, dt);  % user-provided discrete update
        tSol  = [0, dt];
        XSol  = [x0, Xnext];
    else
        %%% <-- CONTINUOUS DYNAMICS (default)
        odeFun = @(t, X) contDynFun(t, X, u);
        [tArr, XArr] = ode45(odeFun, [0, dt], x0);
        tSol = tArr.';
        XSol = XArr.';
    end
end

%% =====================================================================
%% Reconstruct entire (time, state) from z
function [timeVec, XMat] = buildFullTrajectory(z, xInit, uList, contDynFun, discDynFun)
    [tVec, Xbndry] = unpackZ(z, length(uList), length(xInit));
    N = length(uList);

    timeVec = [];
    XMat    = [];

    Xi = xInit;
    t0 = 0;
    for i = 1:N
        dt = tVec(i);
        ui = uList(i);

        [tArc, XArc] = simulateArc(Xi, dt, ui, contDynFun, discDynFun);

        % Shift time
        tArcShifted = t0 + tArc;

        if isempty(timeVec)
            timeVec = tArcShifted;
            XMat    = XArc;
        else
            % avoid repeating last time from previous arc
            timeVec = [timeVec, tArcShifted];
            XMat    = [XMat,   XArc];
        end
        timeVec(end) = timeVec(end)- 1e-8;
        % Next boundary
        Xi = Xbndry(:, i);
        t0 = t0 + dt;
    end
end

%% =====================================================================
%% Piecewise-constant control builder
function uAll = buildControlProfile(timeVec, z, uList, nx)
    [tVec, ~] = unpackZ(z, length(uList), nx);
    tBreaks   = [0; cumsum(tVec)];
    uAll      = zeros(size(timeVec));

    for i = 1:length(uList)
        idx = (timeVec >= tBreaks(i)) & (timeVec < tBreaks(i+1));
        if i == length(uList)
            % Last segment includes the endpoint
            idx = (timeVec >= tBreaks(i)) & (timeVec <= tBreaks(i+1));
        end
        uAll(idx) = uList(i);
    end
end

%% =====================================================================
%% Decision vector packing/unpacking
function z = packZ(tVec, Xbndry)
    % tVec:   (N x 1)
    % Xbndry: (nx x N)
    z = [tVec(:); Xbndry(:)];
end

function [tVec, Xbndry] = unpackZ(z, N, nx)
    % z = [ t1..tN,  X2(1..nx), X3(1..nx), ..., X_{N+1}(1..nx) ]
    tVec = z(1:N);
    xAll = z(N+1:end);
    Xbndry = reshape(xAll, nx, N);
end
