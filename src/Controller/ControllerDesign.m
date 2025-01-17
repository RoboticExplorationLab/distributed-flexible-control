function [OBSWStruct, ModelStruct, IdStruct] = ControllerDesign(OBSWSysDisc, ...
                              IdStruct, ModelStruct, SimParam)
    OBSWStruct       = struct; % Onboard Software Structure
    %% Kalman Decomposition
    % separation of observable/controllable states 
    [TKD, IdKD]    = AdHocKD(OBSWSysDisc(:,IdStruct.IdUWU), ModelStruct, ...
                                                        IdStruct, SimParam);
    OBSWStruct.TKD = TKD;
    IdStruct.IdKD  = IdKD;
    
    %% Guidance Design 
    if ~strcmp(SimParam.GuidMethod,'None')
        NXq      = IdStruct.RealMdl.NXq;
        IdXdrigr = IdStruct.RealMdl.IdXdrigr;
        IdXdmod  = IdStruct.RealMdl.IdXdmod;
        IdXJd    = NXq + [IdXdrigr, IdXdmod];
        MassMat  = ModelStruct.DescMat(IdXJd,IdXJd);
        Mrw      = (ModelStruct.Act.InertiaRW .* ModelStruct.Act.RWDir)';
        OBSWStruct.Guid.MaxTorqueRW = ModelStruct.Act.MaxTorqueRW;
        OBSWStruct.Guid.nomrwspeed = ModelStruct.Act.nomrwspeed;
        OBSWStruct.Guid.maxrwspeed = ModelStruct.Act.maxrwspeed;
        OBSWStruct.Guid.minrwspeed = ModelStruct.Act.minrwspeed;        
        OBSWStruct.Guid.MassMat   = MassMat;
        OBSWStruct.Guid.Mrw       = Mrw;
        OBSWStruct.Guid.InertiaRW = ModelStruct.Act.InertiaRW;
        OBSWStruct.Guid.RWDir     = ModelStruct.Act.RWDir;
    end
    
    %% Feedback Controller Design
    [OBSWStruct, ModelStruct, IdStruct] = LQRDesigner(OBSWSysDisc, ...
                              OBSWStruct, ModelStruct, IdStruct, SimParam);

    %% State Estimation Design
    [OBSWStruct] = PSSKFDesigner(OBSWStruct, OBSWSysDisc, IdStruct,SimParam);
    
end

