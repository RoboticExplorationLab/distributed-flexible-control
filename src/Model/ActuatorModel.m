function [ModelStruct,IdStruct] = ActuatorModel(ModelStruct, IdStruct, SimParam)
    % % 
    SelActId     = ModelStruct.SelActId;
    SelActNodes  = ModelStruct.SelActNodes;
    Node2Mode    = ModelStruct.Node2Mode;
    DampMatModes = ModelStruct.DampMatModes;
    RealId       = IdStruct.RealMdl;

    % % 
    IdStruct.RWperpoint = SimParam.RWperpoint;
    IdStruct.NU         = sum(SimParam.RWperpoint);

    RealId.NXRW   = IdStruct.NU;
    RealId.NXd    = RealId.Nrigd + RealId.Nmod + RealId.NXRW;
    RealId.NX     = RealId.NXq + RealId.NXd;
    RealId.IdXRW  = RealId.Nrig + 2*RealId.Nmod + (1:RealId.NXRW);
    RealId.IdXdRW = RealId.Nrigd + RealId.Nmod + (1:RealId.NXRW);
    RealId.IdXd   = RealId.NXq + (1:RealId.NXd);
    % Each row of RWPos is the position of a Reaction Wheel
    RWPos       = [repelem(SelActNodes,IdStruct.RWperpoint,1)];
    % Each row is the direction of the angular velocity vector of each RW
    RWDir       = SimParam.RWDirPerPoint;
    % Input Jacobian for Linear Model
    BInputJac   = zeros(RealId.NXd,IdStruct.NU);
    % Moment of Inertia of RW [kgm^2]
    InertiaRW   = zeros(IdStruct.NU,1);
    % Max torque [Nm]
    MaxTorqueRW = zeros(IdStruct.NU,1);
    % Dimensions of RW (X8318S) [m] 
    RadRW   = ones(IdStruct.NU,1)*SimParam.RadRW;
    HghtRRW = ones(IdStruct.NU,1)*SimParam.HghtRRW;
    % Input Jacobian link between Reaction Wheel and nodal rotation around
    % rotation axis
    Lrw = zeros(RealId.Nfem,IdStruct.NU);
    % Input Jacobian link between Reaction Wheel and nodal rotation
    Nrw = zeros(IdStruct.NU,3,RealId.Nfem);

    for i=1:IdStruct.NU
        InertiaRW(i)   = SimParam.RWInertia; % 0.1295; % inertia of wheel [kgm^2]
        MaxTorqueRW(i) = SimParam.MaxTorqueRW; % 0.12; % max torque [Nm]
        InertiaRW(i)   = SimParam.RWInertia; % 0.1295; % inertia of wheel [kgm^2]
        MaxTorqueRW(i) = SimParam.MaxTorqueRW; % 0.12; % max torque [Nm]
        % Modal participation matrix of Reaction Wheel
        % Used for the entries of the input jacobian to link reaction wheel
        % actuation to FEM nodal displacement

        % id rot dof
        % get them from tra dof
        % femNrw     = FemN(RWPos(i,:),Nfuncs,ConnectMat,NodeCoord, IdBound);
        % % -<RWDir,[thetax,thetay, 0]>, the latter being the slope of the plate
        % % at the RW position
        if i>SimParam.RWperpoint(1)
            IdPos = find(i<=cumsum(SimParam.RWperpoint),1,'first');
            IdAct = SelActId(IdPos-1)*3 + (-1:0); % (-2:-1);
            a = eye(RealId.Nfem);
            Nrw(i,1:2,:) = a(IdAct,:); 
            Lrw(:,i)   = tensorprod(RWDir(i,:),Nrw(i,:,:),2);
        end
    end 

    Lrw = Node2Mode'*Lrw;
    Nrw = tensorprod(Nrw,Node2Mode,3,1);
    % Input Jacobian: RW effect on nodal displacement and RW speed states
    % DescMatXd\
    BInputJac(RealId.IdXdmod,:) = -Lrw;
    BInputJac(RealId.IdXdRW,:)  = eye(IdStruct.NU);
    BInputJac(RealId.IdXdrigr,:) = -RWDir';

    % StateMat Update
    StateMat = blkdiag(ModelStruct.StateMat, zeros(RealId.NXRW));
    % RW Gyroscopic torque
    Jrwwrw = zeros(3,RealId.Nrigrd + RealId.Nmod);
    Nrwrw  = zeros(RealId.Nmod, RealId.Nrigrd + RealId.Nmod);
    for iu = 1:IdStruct.NU
        jrw0 = skew(RWDir(iu,:)' * (InertiaRW(iu) .* SimParam.nomrwspeed));
        omegarw = squeeze(Nrw(iu,:,:));
        omegarw = [eye(3), omegarw];
        Jrwwrw = Jrwwrw + jrw0 * omegarw;
        Nrwrw  = Nrwrw  + squeeze(Nrw(iu,:,:))' * jrw0 * omegarw;
    end

    StateMat(RealId.IdXrigrd, RealId.IdXrigrd) = Jrwwrw(:,1:3);
    % adding gyricity
    StateMat(RealId.IdXrigrd, RealId.IdXmodd) = Jrwwrw(:,4:end);
    StateMat(RealId.IdXmodd, RealId.IdXrigrd) = Nrwrw(:,1:3);
    StateMat(RealId.IdXmodd, RealId.IdXmodd) = -DampMatModes ...
                                                      + Nrwrw(:,4:end);

    % save params 
    ModelStruct.DescMat            = blkdiag(ModelStruct.DescMat, ...
                                                        diag(InertiaRW));
    ModelStruct.StateMat           = StateMat;
    ModelStruct.Act.MaxTorqueRW    = MaxTorqueRW;
    ModelStruct.Act.RadRW          = RadRW;
    ModelStruct.Act.HghtRRW        = HghtRRW;
    ModelStruct.Act.InertiaRW      = InertiaRW;
    ModelStruct.Act.Lrw            = Lrw;
    ModelStruct.Act.Nrw            = Nrw;
    ModelStruct.Act.BInputJac      = BInputJac;
    ModelStruct.Act.RWDir          = RWDir;
    ModelStruct.Act.RWPos          = RWPos;
    IdStruct.RealMdl               = RealId;
end

