function [ModelStruct,IdStruct] = PerturbationModel(ModelStruct, IdStruct, ...
                                                                SimParam)
    % actuators plus perturbations
    RealId            = IdStruct.RealMdl;
    NpertperRW        = 2*SimParam.RadFPertRW + 2*SimParam.RadMPertRW;
    IdStruct.NW       = NpertperRW * IdStruct.NU;
    IdStruct.NWRW     = NpertperRW * IdStruct.NU;
    IdStruct.NWRWF    = 2*SimParam.RadFPertRW * IdStruct.NU;
    IdStruct.NWRWM    = 2*SimParam.RadMPertRW * IdStruct.NU;
    IdStruct.NUW      = IdStruct.NU + IdStruct.NW;
    IdStruct.IdUWRWU  = 1:IdStruct.NU;
    IdStruct.IdUWRWW  = IdStruct.NU + (1:(NpertperRW * IdStruct.NU));
    IdStruct.IdUWRWFW = IdStruct.NU + (0:NpertperRW:(NpertperRW * (IdStruct.NU-1))) ...
                                                + (1:2*SimParam.RadFPertRW)';
    IdStruct.IdUWRWFW = IdStruct.IdUWRWFW(:)';
    IdStruct.IdUWRWMW = IdStruct.NU + (0:NpertperRW:(NpertperRW * (IdStruct.NU-1))) ...
                        + SimParam.RadFPertRW*2 + (1:2*SimParam.RadMPertRW)'; 
    IdStruct.IdUWRWMW = IdStruct.IdUWRWMW(:)';
    IdStruct.IdUWW    = IdStruct.NU + (1:IdStruct.NW);
    IdStruct.IdUWWRW  = 1:(NpertperRW * IdStruct.NU);
    IdStruct.IdUWWRWF = IdStruct.IdUWRWFW - IdStruct.NU;
    IdStruct.IdUWWRWM = IdStruct.IdUWRWMW - IdStruct.NU;
    IdStruct.IdUWU    = 1:IdStruct.NU;
    IdStruct.IdUWURW  = 1:IdStruct.NU;

    SelActId    = ModelStruct.SelActId;
    Node2Mode        = ModelStruct.Node2Mode;
    GDistJac         = zeros(RealId.NX,IdStruct.NW);
    ModelStruct.Pert = struct;
    
    % if SimParam.rigidmodes
    IdBound    = ModelStruct.IdBound*3 + (-2:0);
    IdBound    = IdBound(:);
    ConnectMat = ModelStruct.ConnectMat;
    NodeCoord  = ModelStruct.NodeCoord;
    Nfuncs     = ModelStruct.Nfuncs;

    % % Jitter perturbations
    % RW Jitter
    % broadband noise level x2 (Force and Torque)
    % First Harmonic Coefficient x2 (F and M)
    if SimParam.RadFPertRW
        ModelStruct.Pert.FRWRadBBVar   = SimParam.FRWRadBBVar;
        ModelStruct.Pert.FRWRadFHCoeff = SimParam.FRWRadFHCoeff;
    else
        ModelStruct.Pert.FRWRadBBVar   = 0;
        ModelStruct.Pert.FRWRadFHCoeff = 0;
    end
    if SimParam.RadMPertRW
        ModelStruct.Pert.MRWRadBBVar   = SimParam.MRWRadBBVar;
        ModelStruct.Pert.MRWRadFHCoeff = SimParam.MRWRadFHCoeff;
    else
        ModelStruct.Pert.MRWRadBBVar   = 0;
        ModelStruct.Pert.MRWRadFHCoeff = 0;
    end

    RWDir       = ModelStruct.Act.RWDir;
    RWPos       = ModelStruct.Act.RWPos;
    Nrw         = ModelStruct.Act.Nrw;
    IdXNLdw     = [RealId.IdXrigd, RealId.IdXmodd];
    if NpertperRW
        for iu=1:IdStruct.NU
            IdWFXY   = NpertperRW*(iu-1) + (1:2);
            IdWMXY   = NpertperRW*(iu-1) + (3:4);
            RWRadDir = null(RWDir(iu,:));
            
            if iu>SimParam.RWperpoint(1)
                IdPos = find(iu<=cumsum(SimParam.RWperpoint),1,'first');
                IdAct = SelActId(IdPos-1)*3 + -2; % (-2:-1);
                Etrw  = Node2Mode(IdAct,:)'*RWRadDir(3,:);
            else
                Etrw = zeros(IdStruct.RealMdl.Nmod,2);
            end
            if SimParam.RadFPertRW
                % X,Y Radial Jitter Force ModelStruct.Pert.FRWRadFHCoeff
                GDistJac(IdXNLdw,IdWFXY) =               [RWRadDir; % rigt
                                    skew([RWPos(iu,:),0])*RWRadDir; % rigr
                                                             Etrw]; % mod
            end
            if SimParam.RadMPertRW
                % X,Y Radial Jitter Momentum
                GDistJac(RealId.IdXrigrd,IdWMXY) = RWRadDir;
                % ModelStruct.Pert.MRWRadFHCoeff * RWRadDir;
            end
    
            Errw     = squeeze(Nrw(iu,:,:))'*RWRadDir;
            GDistJac(RealId.IdXmodd,IdWMXY) = Errw;
        end
    end

    WPert = eye(IdStruct.NWRW);
    WPert(IdStruct.IdUWWRWF,:) = ModelStruct.Pert.FRWRadFHCoeff*WPert(IdStruct.IdUWWRWF,:);
    WPert(IdStruct.IdUWWRWM,:) = ModelStruct.Pert.MRWRadFHCoeff*WPert(IdStruct.IdUWWRWM,:);
    
    WPertUnkown = eye(IdStruct.NWRW);
    WPertUnkown(IdStruct.IdUWWRWF,:) = ModelStruct.Pert.FRWRadBBVar*WPertUnkown(IdStruct.IdUWWRWF,:);
    WPertUnkown(IdStruct.IdUWWRWM,:) = ModelStruct.Pert.MRWRadBBVar*WPertUnkown(IdStruct.IdUWWRWM,:);
    ModelStruct.Pert.WPert    = WPert;
    ModelStruct.Pert.WPertUnkown = WPertUnkown;
    ModelStruct.Pert.GDistJac = GDistJac;

end

