function [OBSWStruct] = PSSKFDesigner(OBSWStruct, SysIddModel, ...
                                                         IdStruct,SimParam)
    %% Load Vars
    IdKD = IdStruct.IdKD;
    IdXRW = [IdKD.IdXNullRWd, IdKD.IdXMomRW];
    
    %% Feedback Design
    % % disregarding wheel speeds
    Nxx     = IdStruct.RealMdl.NX;
    Nyy     = IdStruct.NY;
    Q = SimParam.QEst;
    R = SimParam.REst;

    % % define an observer only for detectable subspace
    % % z = T * x; x = T'z;
    nxest = IdKD.NOB;
    T2Est = OBSWStruct.TKD(:,IdKD.IdOB)';
    TT = (T2Est' * T2Est);
    Ad = T2Est * SysIddModel.A * T2Est';
    Cd = SysIddModel.C * T2Est';
    Dd = SysIddModel.D(:,IdStruct.IdUWU);
    
    NT   = 10000;
    Lall = zeros(Nxx, Nyy, SimParam.NStepSST);
    P    = zeros(nxest, nxest);
    ptrace = zeros(1, NT);
    for i=2:NT
        ikf = 1+mod(i-1,SimParam.NStepSST);
        if ikf == 1
            IdY = 1:IdStruct.NY;
        else
            IdY = setdiff(1:IdStruct.NY, IdStruct.IdYSST);
        end
        Ri  = R(IdY,IdY);
        Ci  = Cd(IdY,:);
        P = Ad * P * Ad' + Q;
        % Update step in Joseph's form
        S = Ci * P * Ci' + Ri;
        K = P * Ci' / S;
        M = (eye(size(Ad)) - K * Ci);
        P      = M * P * M' + K * Ri * K';
        Lall(:,IdY,ikf) = T2Est' * K;
        ptrace(i) = trace(P);
    end


    %% Save vars
    OBSWStruct.Est        = struct;
    OBSWStruct.Est.NxEst  = 1:nxest;
    OBSWStruct.Est.IdXRW  = IdXRW;
    OBSWStruct.Est.T2Est  = T2Est;
    OBSWStruct.Est.IdXEst = IdKD.IdOB;
    OBSWStruct.Est.L      = Lall; 
    OBSWStruct.Est.A      = T2Est' * Ad * T2Est;
    OBSWStruct.Est.B      = TT * SysIddModel.B(:,IdStruct.IdUWU);
    OBSWStruct.Est.C      = Cd * T2Est;
    OBSWStruct.Est.D      = Dd;
end

