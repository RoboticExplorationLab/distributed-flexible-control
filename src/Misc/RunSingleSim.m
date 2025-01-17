function ResStruct = RunSingleSim(ModelStruct,OBSWStruct, RealStruct, ...
                                                        IdStruct, SimParam)
    SamplingTime = ModelStruct.Sens.SamplingTime;
    IdXRW        = IdStruct.RealMdl.IdXRW;
    RealID       = IdStruct.RealMdl;
    NX           = RealID.NX;
    NU           = IdStruct.NU;
    NW           = IdStruct.NW;
    nom_rw       = ModelStruct.Act.nomrwspeed;
    tVec         = 0:SamplingTime:SimParam.tSim;

    NT = length(tVec);
    % full system matrices
    % xk+1=AAxk + BBuk + GGwk
    % yk = CCxk + DDuk
    RealSysDisc = RealStruct.RealSysDisc;
    GG = RealSysDisc.B(:,IdStruct.IdUWW);

    % % states
    % closed loop states
    CLStates      = zeros(NX,NT);
    % open-loop states
    OLStates      = CLStates;
    % initial state
    x_0               = zeros(NX,1);
    x_0(RealID.IdXRW) = nom_rw; 
    x_0(RealID.IdXrigrq) = (SimParam.InitAttErr.^0.5)/4;
    CLStates(:,1)     = x_0;
    OLStates(:,1)     = x_0;
    x_nom             = zeros(NX,1);
    x_nom(IdXRW)   = nom_rw;

    thetarwCL = zeros(RealID.NXRW,1);
    thetarwOL = thetarwCL;
    
    % % measurements
    % closed-loop real world model measurements
    CLMeas = zeros(IdStruct.NY,NT);
    % open-loop real world model measurements
    OLMeas = zeros(IdStruct.NY,NT);
    
    UVec  = zeros(NT,NU);
    x_Est = zeros(NX,NT);
    x_Est(IdXRW,1) = nom_rw;
    y_hat = zeros(IdStruct.NY,NT);
    wKnown   = zeros(NW,NT);
    wUnknown = zeros(NW,NT);
    wPert = zeros(NW,NT);
    % Slewing Guidance Profile
    [x_ref, u_ref] = StateGuidance( tVec, RealStruct.RealSysCont,...
                                        SimParam, OBSWStruct, IdStruct);

    for it = 1:NT-1
        if it == 1
            u_prev     = zeros(NU,1);
            w_known_prev = zeros(NW,1);
            x_hat_prev = x_Est(:,1);
        else
            u_prev     = UVec(it-1,:)';
            w_known_prev = wKnown(:,it-1);
            x_hat_prev = x_Est(:,it-1);
        end

        % perturbation: jitter
        % split rw disturbances between forces and moments
        if SimParam.SimRWPertFlag
            [wKnown(:,it),wUnknown(:,it)] = RWPertFunc(thetarwCL,CLStates(IdXRW,it),...
                                            ModelStruct,IdStruct,SimParam); 
            wPert(:,it) = wKnown(:,it) + wUnknown(:,it);
            [wKnownOL,wUnknownOL] = RWPertFunc(thetarwOL,OLStates(IdXRW,it),...
                                            ModelStruct,IdStruct,SimParam);
            wPertOL = (wKnownOL + wUnknownOL)';
        else
            wPertOL = zeros(NW,1);
        end
        
        
        % measurement function 
        [CLMeas(:,it),~] = ObsFunc(RealStruct, CLStates(:,it), u_prev, ...
                    wPert(:,it), thetarwCL, IdStruct, ModelStruct, SimParam);
        [OLMeas(:,it),~] = ObsFunc(RealStruct, OLStates(:,it), zeros(NU,1), ...
                        wPertOL, thetarwOL, IdStruct, ModelStruct, SimParam);

        % State Estimation (State Feedback)
        [x_Est(1:NX,it), y_hat(:,it)] = ...
                StateEstimation(OBSWStruct.Est, CLMeas(:,it), u_prev, w_known_prev, ...
                x_nom, x_hat_prev, GG, SimParam, it);
        
        % Update initial state 
        UVec(it,:) = StateControl(x_nom, x_ref, u_ref, x_Est(:,it), it, ...
                                     OBSWStruct.Ctrl, IdStruct, SimParam);

        disp(['elapsed time: ' num2str(tVec(it)) ' s/total: ' num2str(tVec(end)) ' s'])
            
        % propagation function
        [CLStates(:,it+1),thetarwCL] = ...
            rigflexdynDisc(CLStates(:,it), UVec(it,:)', thetarwCL, ...
                                RealStruct, IdStruct, ModelStruct, SimParam);

        [OLStates(:,it+1),thetarwOL] =...
            rigflexdynDisc(OLStates(:,it),zeros(NU,1), thetarwOL, ...
                                RealStruct, IdStruct, ModelStruct, SimParam);
    end
    
    %% save results to output
    ResStruct = struct;
    % add reference values again (only applicable to linear system)
    CLStates = CLStates';
    OLStates = OLStates';
    ResStruct.tVec       = tVec;
    ResStruct.y_hat      = y_hat;
    ResStruct.UVec       = UVec;
    ResStruct.x_hat      = x_Est;
    ResStruct.CLMeas     = CLMeas;
    ResStruct.CLStates   = CLStates;
    ResStruct.OLStates   = OLStates;
    ResStruct.x_ref      = x_ref;
    ResStruct.u_ref      = u_ref;
    ResStruct.RealStruct = RealStruct;
end

