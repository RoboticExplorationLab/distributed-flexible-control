function [x_ref, u_ref] = ZVTORigSlewGuidance( tVec, RealSysCont,...
                                        SimParam, OBSWStruct, IdStruct)
    % load vars
    t_start     = SimParam.TSlewStart;
    Att_Tgt     = zeros(3,1);
    SlewRng     = SimParam.TgtSlewAngleRng;

    Att_Tgt(1)   = SlewRng(2);
    IdXrigrq     = IdStruct.RealMdl.IdXrigrq(1);
    IdXrigrd     = IdStruct.RealMdl.IdXrigrd(1);

    GuidStruct   = OBSWStruct.Guid;
    MassMat      = GuidStruct.MassMat;
    Mrw          = GuidStruct.Mrw;
    InertiaRW    = OBSWStruct.Guid.InertiaRW;
    nomrwspeed   = GuidStruct.nomrwspeed;
    maxrwspeed   = GuidStruct.maxrwspeed;
    minrwspeed   = GuidStruct.minrwspeed;
    MaxTorqueRW  = GuidStruct.MaxTorqueRW;
    NT           = length(tVec);
    NX           = IdStruct.RealMdl.NX;
    IdKD         = IdStruct.IdKD;
    NU           = IdStruct.NU;
    % durations
    Jsc   = MassMat(1:3,1:3);
    HCap  = abs(Mrw) * (maxrwspeed-nomrwspeed);
    Hslew = HCap*(4/5);

    IdX = [IdXrigrq, IdXrigrd, IdKD.IdXMomRW(1)];
    TOKD = OBSWStruct.TKD(:,IdX);
    TOKD = TOKD./sum(abs(TOKD),1);
    Tact  = OBSWStruct.TKD(IdStruct.RealMdl.IdXRW,[IdKD.IdXMomRW])';
    Tact = Tact./vecnorm(Tact,2,1);
    umax = abs(Tact) * MaxTorqueRW*4/5;

    % Build the data struct
    problemStruct.xInit   = zeros(3,1);     % initial state
    problemStruct.xTarget = [Att_Tgt(1)*pi/180;0;0];   % final state
    problemStruct.uList   = [-umax(1);0;umax(1)];     % bang-coast-bang controls
    problemStruct.guessT = [2.5, 95, 2.5];     % initial guess for durations
    problemStruct.guessX = zeros(3,3); % guess boundary states X2..X4
    problemStruct.dynFun = @(t, X, u) singleAxisDynamics(X, u, Jsc(1,1));
    problemStruct.arcConstraintFun = @(tArray, XArray) momentumConstraint(XArray, Hslew(1));
    problemStruct.costFun = @(tVec, XMat) tVec(end);
    problemStruct.multistart = false;
    problemStruct.NMStart    = 10;
    problemStruct.maxT       = SimParam.tSim;

    % optimize bang-coast-bang profile
    [timeVec, XVec, UVec] = solveMultiArcManeuverMS(problemStruct, []);

    % ZVD filtering
    % get the natural frequencies and poles to filter out
    NewB = RealSysCont.B(:,1:NU) * pinv(Tact);
    newC = eye(NX);  newC = newC(IdXrigrq(1),:);
    newss = ss(RealSysCont.E\RealSysCont.A, RealSysCont.E\NewB(:,1), newC,0);
    [msys,U] = minreal(newss,eps*100); % sqrt(eps)*100);
    % sysT = ss2ss(sys,T)
    U = U(1:length(msys.StateName),:);
    R = reducespec(msys,"balanced");
    R.Options.Goal = "relative";
    R = process(R);
    t = 5e-2; % 2; % 5e-2; % 
    Nx = find(R.Energy<t,1,'first')-1;
    [rsys,info] = getrom(R,Order=Nx,Method="truncate");
    sumnum = zeros(1,5*SimParam.SamplingFreq);
    for i=1:(5*SimParam.SamplingFreq)
        [wn,zeta,P] = damp(c2d(rsys,i/SimParam.SamplingFreq));
        P = P(abs(abs(P)-1) > 1e-5);
        ZVD = zpk(P,zeros(1,length(P)),1,1/SimParam.SamplingFreq);
        [num,den] = tfdata(ZVD,'v');
        sumnum(i) = min(num);
        if sumnum(i) > 0
            ZVD =  ZVD / sum(num);
            [num,den] = tfdata(ZVD,'v');
            filtzvd = zeros(i*(length(P)+1),1);
            filtzvd(1:i:end) = num;
            break;
        end
    end
    % post-process results
    % interpolate timeVec to tVec
    % need index relating states in XVec to x_ref
    timeVec = timeVec' + t_start;
    timeVec = [0; t_start-1e2*eps; timeVec; timeVec(end)+1e2*eps; 1e10];
    XVec    = [XVec(:,1), XVec(:,1), XVec, XVec(:,end), XVec(:,end)];
    UVec    = [0 0 UVec 0 0];

    x_ref   = zeros(NX, NT);
    u_ref   = zeros(NU, NT);
    InvTact = pinv(Tact);
    for i=1:size(UVec,1)
        UVecnewT = interp1(timeVec,UVec(i,:),tVec);
        UVecnewTFilt = conv(UVecnewT,filtzvd ,'same');
        u_ref = u_ref + InvTact(:,i)*UVecnewTFilt;
    end
    % roll-forward
    IdXrigrq = IdStruct.RealMdl.IdXrigrq;
    IdXrigrd = IdStruct.RealMdl.IdXrigrd;
    IdXmodd   = IdStruct.RealMdl.IdXmodd;

    if SimParam.etaref
        A = RealSysCont.A;
        A(IdXrigrd, [IdXrigrd, IdXmodd]) = 0;
        A([IdXrigrd, IdXmodd], IdXrigrd) = 0;
        A(IdXmodd, IdXmodd) = diag(diag(A(IdXmodd, IdXmodd)));
        RealSysCont.A = A;
        RealSysDisc = c2d(RealSysCont,1/SimParam.SamplingFreq,'zoh');
        Ad = RealSysDisc.A; Bd = RealSysDisc.B(:,1:NU);
        for i=2:NT
            x_ref(:,i) = Ad*x_ref(:,i-1) + Bd*u_ref(:,i-1);
        end
    else
        XVecnewT = interp1(timeVec,XVec(1,:),tVec);
        x_ref(IdXrigrq(1),:) = XVecnewT/4;
        XVecnewT = interp1(timeVec,XVec(2,:),tVec);
        x_ref(IdXrigrd(1),:) = XVecnewT;
        XVecnewT = interp1(timeVec,XVec(3,:),tVec);
        x_ref = x_ref + TOKD(:,3)*XVecnewT;
        x_ref(IdStruct.RealMdl.IdXRW,:) = x_ref(IdStruct.RealMdl.IdXRW,:)./InertiaRW;
    end
    %{%}

    % to correct numerical precision attitude drift 
    % in guidance profile (order of 0.0002 deg/0.7" for 0.5 deg slew)
    cor_factor = x_ref(IdXrigrq(1),end)*4*180/pi/Att_Tgt(1);
    x_ref(IdXrigrq(1),:) = x_ref(IdXrigrq(1),:) / cor_factor;

    close all;
    IdXdrigrq = IdStruct.RealMdl.IdXrigrq;
    ylbls = {'X [deg]','Y [deg]','Z [deg]'};
    figure(1);
    for i=1:3
        subplot(3,1,i)
        plot(tVec,x_ref(IdXdrigrq(i),:)*4*180/pi);
        grid on;
        ylabel(ylbls{i})
        if i ==1
            title('SC Rotation Angle')
        end
    end

    IdXdrigrd = IdStruct.RealMdl.IdXrigrd;
    figure(2);
    ylbls = {'X [rad/s]','Y [rad/s]','Z [rad/s]'};
    for i=1:3
        subplot(3,1,i)
        plot(tVec,x_ref(IdXdrigrd(i),:)*180/pi);
        grid on;
        ylabel(ylbls{i},'fontsize',12)
        if i==1
            title('SC Angular velocity')
        end
    end
    xlabel('time [s]','fontsize',12)

    Nrwpp = min(3,max(IdStruct.RWperpoint));
    NT = length(tVec);
    figure(3);
    ylbls = {'X [rpm]','Y [rpm]','Z [rpm]'};
    for j=1:Nrwpp
        subplot(Nrwpp,1,j)
        subplot(Nrwpp,1,j)
        hold on;
        plot(tVec,(maxrwspeed-nomrwspeed)*60/(2*pi)*ones(1,NT),'k--');
        plot(tVec,(minrwspeed-nomrwspeed)*60/(2*pi)*ones(1,NT),'k--');
        
        for i=1:length(IdStruct.RWperpoint)
            Nrwm1 = sum(IdStruct.RWperpoint(1:i-1));
            if IdStruct.RWperpoint(i) >= j
                wrw = x_ref(IdStruct.RealMdl.IdXRW(Nrwm1+j),:)'*60/(2*pi);
                dwrw = wrw; %  - nomrwspeed(Nrwm1+j)*60/(2*pi);
                plot(tVec,dwrw,'k')
            end
        end
        grid on;
        ylabel(ylbls{j},'fontsize',12)
        if j ==1
            title('RW Speed wrt nominal Speed')
        end
    end
    xlabel('time [s]','fontsize',12)

    figure(4);
    ylab = {'X [Nm]','Y [Nm]','Z [Nm]'};
    % assuming all RW are facing x,y,z at every node in order
    for j= 1:min(3,max(IdStruct.RWperpoint))
        subplot(min(3,max(IdStruct.RWperpoint)),1,j)
        hold on;
        plot(tVec,MaxTorqueRW(1)*ones(1,NT),'k--');
        plot(tVec,-MaxTorqueRW(1)*ones(1,NT),'k--');
        for i=1:length(IdStruct.RWperpoint)
            if IdStruct.RWperpoint(i) >= j
                stairs(tVec,u_ref(sum(IdStruct.RWperpoint(1:i-1))+j,:)','k')
            end
        end
        grid on;
        ylabel(ylab{j})
        if j == 1
            title('RW torque Commands')
        end
    end
    xlabel('time [s]')
end

function dXdt = singleAxisDynamics(X, u, J)
    % Single-axis with reaction wheel, as an example
    theta = X(1);
    omega = X(2);
    h     = X(3);

    dtheta = omega;
    domega = -(1/J)*u;
    dh     = u;
    dXdt   = [dtheta; domega; dh];
end

function [c, ceq] = momentumConstraint(XArray, hmax)
    % Example arc constraint: ensure |h(t)| <= hmax over the arc
    hAll = XArray(3,:);
    c    = max(abs(hAll)) - hmax;  % must be <= 0
    ceq  = [];
end