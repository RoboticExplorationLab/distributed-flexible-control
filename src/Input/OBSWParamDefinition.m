function SimParam = OBSWParamDefinition(SimParam, SysIddModel, ModelStruct, ...
                                        IdStruct, TKD, IdKD, ActPlacementFlag, ...
                                        ScenarioFlag, GuidProfileFlag)

    if strcmp(ActPlacementFlag,'Centralized')
        ScaleFactor = 11/2;
    elseif strcmp(ActPlacementFlag,'Distributed')
        ScaleFactor = 1;
    else
        msg = ['Actuator Placement Flag ' ActPlacementFlag ' not valid.' ...
                ' Valid Options are "Centralized" and "Distributed'];
        error(msg)
    end

    if strcmp(ScenarioFlag,'Slewing')
        SimParam.GuidMethod = 'ZVTOrig'; 
    elseif strcmp(ScenarioFlag,'FinePointing')
        SimParam.GuidMethod = 'None'; 
    else
        msg = ['Scenario Flag ' ScenarioFlag ' not valid.' ...
                ' Valid Options are "Slewing" and "FinePointing'];
        error(msg)
    end

    if strcmp(GuidProfileFlag,'TimeOptimal')
        SimParam.etaref     = true;
    elseif strcmp(GuidProfileFlag,'Smoothed')
        SimParam.etaref     = false;
    else
        msg = ['Guidance Profile Flag ' GuidProfileFlag ' not valid.' ...
                ' Valid Options are "TimeOptimal" and "Smoothed'];
        error(msg)
    end
    SimParam.Nhor     = 100;

    % % Controller tuning 
    NxCtrl = IdKD.NCO;
    NuCtrl = IdStruct.NU;
    
    Q = zeros(NxCtrl); 
    Q(IdKD.IdXFiltHPAMU,IdKD.IdXFiltHPAMU) = 1e4*eye(IdKD.NFiltHPAMU);
    Q(IdKD.IdXFiltIMU,IdKD.IdXFiltIMU) = 1e4*eye(IdKD.NFiltIMU);
    Q(IdKD.IdXFiltSST,IdKD.IdXFiltSST) = 1e10*eye(IdKD.NFiltSST);
    Q(IdKD.IdXrigrq,IdKD.IdXrigrq)     = 1e7*eye(IdKD.Nrigrq);
    Q(IdKD.IdXmodq,IdKD.IdXmodq)       = 1e8*ModelStruct.StiffMatModes^(-1);
    Q(IdKD.IdXrigrd,IdKD.IdXrigrd)     = 1e3*eye(IdKD.Nrigrd);
    Q(IdKD.IdXmodd,IdKD.IdXmodd)       = 1e1*ModelStruct.StiffMatModes^(-1);
    WRW = 1e-1/ScaleFactor; 
    Q(IdKD.IdXNullRWd,IdKD.IdXNullRWd) = WRW*eye(IdKD.NNullRWd); 
    
    R = 1e-1*eye(NuCtrl)/ScaleFactor;
    
    QN = 1*Q;
    
    SimParam.QCtrl  = Q;
    SimParam.RCtrl  = R;
    SimParam.QNCtrl = QN;
    
    % distributed
    NxEst = IdKD.NOB;
    % tradeoff between robustness to uncertainty and error convergence/speed
    Gd = TKD(:,IdKD.IdOB)' * SysIddModel.B(:,IdStruct.IdUWWRW);
    QEst  = 1e-12*eye(NxEst) + Gd * ModelStruct.Pert.WPertUnkown * Gd';
    
    REst = diag([SimParam.GyroWNVar*ones(1,IdStruct.NYGyro),...
                  SimParam.AMUWNVar*ones(1,IdStruct.NYAMU),...
                  SimParam.RWWNVar*ones(1,IdStruct.NYRW),...
                  SimParam.SSTWNVar*ones(1,IdStruct.NYSST)]); 
    SimParam.QEst = QEst;
    SimParam.REst = REst;
end

