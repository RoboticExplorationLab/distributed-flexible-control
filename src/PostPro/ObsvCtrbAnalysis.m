function [TKD, IdKD] = ObsvCtrbAnalysis(OBSWSysDisc,ModelStruct, ...
                                                       IdStruct, SimParam)
    OBSWSysDiscnoRW = OBSWSysDisc(:,IdStruct.IdUWU);
    OBSWSysDiscnoRW = xelim(OBSWSysDiscnoRW,IdStruct.RealMdl.IdXRW,"Truncate");
    OBSWSysDisc0 = OBSWSysDisc(:,IdStruct.IdUWU);
    [GS,GNS] = stabsep(OBSWSysDisc0);
    St = size(GNS.a,1) - rank(ctrb(GNS.a,GNS.b));
    disp(['unstabilizable states: ' num2str(St)])
    % unstabilizable states should be number of rigid body modes
    % would need adding actuators such as thrusters too make everything
    % stabilizable
    Co = size(OBSWSysDisc0.A,1) - rank(ctrb(OBSWSysDisc0.A,OBSWSysDisc0.B));
    disp(['uncontrollable states: ' num2str(Co)])
    Ob = size(OBSWSysDisc0.A,1) - rank(obsv(OBSWSysDisc0.A,OBSWSysDisc0.C));
    disp(['unobservable states: ' num2str(Ob)])
    % notwithstanding the rw speed, the system is stabilizable and observable
    
    disp('without rw speed: ')
    [GS,GNS] = stabsep(OBSWSysDiscnoRW);
    St = size(GNS.a,1) - rank(ctrb(GNS.a,GNS.b));
    disp(['unstabilizable states: ' num2str(St)])
    % unstabilizable states should be number of rigid body modes
    % would need adding actuators such as thrusters too make everything
    % stabilizable
    Co = size(OBSWSysDiscnoRW.A,1) - rank(ctrb(OBSWSysDiscnoRW.A,OBSWSysDiscnoRW.B));
    disp(['uncontrollable states: ' num2str(Co)])
    Ob = size(OBSWSysDiscnoRW.A,1) - rank(obsv(OBSWSysDiscnoRW.A,OBSWSysDiscnoRW.C));
    disp(['unobservable states: ' num2str(Ob)])
    % notwithstanding the rw speed, the system is stabilizable and observable
    
    disp('without rw non-nullspace speed and linear velocity:')
    % change of state space of actuation
    % % define a controller only for reachable states
    [TKD, IdKD]    = AdHocKD(OBSWSysDisc(:,IdStruct.IdUWU), ModelStruct, ...
                                                            IdStruct, SimParam);
    OBSWSysDiscCtrOnly = ss2ss(OBSWSysDisc(:,IdStruct.IdUWU),TKD');
    RWCtrMat = TKD(IdStruct.RealMdl.IdXRW, IdKD.IdXRW);
    OBSWSysDiscCtrOnly = OBSWSysDiscCtrOnly * RWCtrMat;
    OBSWSysDiscCtrOnly = xelim(OBSWSysDiscCtrOnly,IdKD.IdUC,"truncate");
    [GS,GNS] = stabsep(OBSWSysDiscCtrOnly);
    St = size(GNS.a,1) - rank(ctrb(GNS.a,GNS.b));
    disp(['unstabilizable states: ' num2str(St)])
    % unstabilizable states should be number of rigid body modes
    % would need adding actuators such as thrusters too make everything
    % stabilizable
    Co = size(OBSWSysDiscCtrOnly.A,1) - rank(ctrb(OBSWSysDiscCtrOnly.A,OBSWSysDiscCtrOnly.B));
    disp(['uncontrollable states: ' num2str(Co)])
    Ob = size(OBSWSysDiscCtrOnly.A,1) - rank(obsv(OBSWSysDiscCtrOnly.A,OBSWSysDiscCtrOnly.C));
    disp(['unobservable states: ' num2str(Ob)])
    
    % % define a controller only for reachable states
    % % z = T * x; x = T'z;
end

