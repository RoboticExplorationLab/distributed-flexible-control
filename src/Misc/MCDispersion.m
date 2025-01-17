function RealStruct = MCDispersion(ModelStruct, RealSysCont, IdStruct, SimParam)
    %MCDISPERSION Summary of this function goes here
    %   Detailed explanation goes here
    RealStruct = struct;
    RealID     = IdStruct.RealMdl;
    Nmod   = RealID.Nmod;
    Nrigd  = RealID.Nrigd;
    Nrigtd = RealID.Nrigtd;
    Nrigrd = RealID.Nrigrd;
    NXRW   = RealID.NXRW;
    Nrw    = ModelStruct.Act.Nrw;
    NFilt  = IdStruct.RealMdl.NFilt;

    % Disperse bias
    ModelStruct.Sens.GyroBias = SimParam.GyroBiasVar*randn(IdStruct.NYGyro,1);
    ModelStruct.Sens.AMUBias  = SimParam.AMUBiasVar*randn(IdStruct.NYAMU,1);
    
    % Disperse parameters
    InertiaCov       = SimParam.InertiaCov;
    ModPartRotMatCov = SimParam.ModPartRotMatCov;
    ModPartLinMatCov = SimParam.ModPartLinMatCov;
    omegaCov         = SimParam.omegaCov;
    dampCov          = SimParam.dampCov;
    BInputJacCov     = SimParam.BInputJacCov;
    MeasMatCov       = SimParam.MeasMatCov;
    
    J  = ModelStruct.SCInertia*diag(1+InertiaCov.^0.5*randn(1,3));
    m  = ModelStruct.SCMass;
    Lr = ModelStruct.ModPartRotMat.*(1+...
        ModPartRotMatCov.^0.5*randn(size(ModelStruct.ModPartRotMat)));
    Lt = ModelStruct.ModPartTransMat.*(1+...
        ModPartLinMatCov.^0.5*randn(size(ModelStruct.ModPartTransMat)));
    Jrw   = ModelStruct.Act.InertiaRW;
    RWDir = ModelStruct.Act.RWDir;
    omega = diag(ModelStruct.StiffMatModes).^0.5;
    xeta  = diag(ModelStruct.DampMatModes)./(2*omega);
    omega = diag(1+omegaCov.^0.5*randn(1,Nmod))*omega;
    K     = diag(omega.^2);
    xeta  = diag(1+dampCov.^0.5*randn(1,Nmod))*xeta;
    C     = diag(2.*xeta.*omega);

    % Perturbation jacobian (not considering model perturbation).
    G = ModelStruct.Pert.GDistJac(RealID.IdXd-NFilt,IdStruct.IdUWWRW);
    
    M = zeros(Nrigd+Nmod+NXRW);
    M(1:3,1:3) = m*eye(Nrigtd); M(4:6,4:6) = J;
    M(1:3,6+(1:Nmod)) = Lt; M(4:6,6+(1:Nmod)) = Lr;
    M(6+(1:Nmod),6+(1:Nmod)) = eye(Nmod);
    M(6+Nmod + (1:NXRW),6+Nmod + (1:NXRW)) = diag(Jrw);
    M(6+(1:Nmod),1:3) = Lt';
    M(6+(1:Nmod),4:6) = Lr';

    Minv = inv(M);
    B = ModelStruct.Act.BInputJac;
    B(RealID.IdXdmod,:) = B(RealID.IdXdmod,:).*(1+BInputJacCov.^0.5*randn(Nmod,NXRW));
    
    RealStruct.m     = m;
    RealStruct.J     = J;
    RealStruct.Leta  = Lr;
    RealStruct.Lr    = Lt;
    RealStruct.Jrw   = Jrw;
    RealStruct.Nrw   = Nrw;
    RealStruct.EDir  = ModelStruct.Act.RWDir;
    RealStruct.Minv  = Minv;
    RealStruct.MinvB = Minv*B;
    RealStruct.MinvK = Minv(:,RealID.IdXdmod)*K;
    RealStruct.MinvC = Minv(:,RealID.IdXdmod)*C;
    RealStruct.MinvG = Minv*G;
    RealStruct.SampTime = 1/SimParam.SamplingFreq;

    % State Mat
    IdXmodq = RealID.IdXmodq;
    IdXmodd = RealID.IdXmodd;
    StateMat = zeros(RealID.NX);
    StateMat(IdXmodq, IdXmodd) = eye(Nmod);
    StateMat(IdXmodd, IdXmodq) = -K;
    StateMat(IdXmodd, IdXmodd) = -C;
    % integration of rigid body motion
    IdXrigtq = RealID.IdXrigtq;
    IdXrigrq = RealID.IdXrigrq;
    IdXrigtd = RealID.IdXrigtd;
    IdXrigrd = RealID.IdXrigrd;
    StateMat(IdXrigtq, IdXrigtd) = eye(3);
    StateMat(IdXrigrq, IdXrigrd) = eye(3)/4;
    % RW Gyroscopic torque
    IdXrigrd = RealID.IdXrigrd;
    IdXFilt  = RealID.IdXFilt;
    Nrigrd   = RealID.Nrigrd;
    Nmod     = RealID.Nmod;
    NX       = RealID.NX;
    % StateMat(IdXrigrd, IdXrigrd) = skew(RWDir' * (Jrw .* ModelStruct.Act.nomrwspeed));

    Jrwwrw = zeros(3,Nrigrd + Nmod);
    Nrwrw  = zeros(Nmod, Nrigrd + Nmod);
    for iu = 1:IdStruct.NU
        jrw0    = skew(RWDir(iu,:)' * (Jrw(iu) .* ModelStruct.Act.nomrwspeed(iu)));
        omegarw = [eye(3), squeeze(Nrw(iu,:,:))];
        Jrwwrw  = Jrwwrw + jrw0 * omegarw;
        Nrwrw   = Nrwrw  + squeeze(Nrw(iu,:,:))' *  jrw0 * omegarw;
    end
    % adding gyricity to linear model
    StateMat(IdXrigrd, IdXrigrd) = Jrwwrw(:,1:3);
    StateMat(IdXrigrd, IdXmodd)  = Jrwwrw(:,4:end);
    StateMat(IdXmodd, IdXrigrd)  = Nrwrw(:,1:3);
    StateMat(IdXmodd, IdXmodd)   = -C + Nrwrw(:,4:end);
    StateMat(IdXFilt,:)         = [];
    StateMat(:,IdXFilt)         = [];

    % Input Mat
    InputMat                             = zeros(NX, IdStruct.NUW);
    InputMat(RealID.IdXd,IdStruct.IdUWU) = B; 
    InputMat(:,IdStruct.IdUWRWW)         = RealSysCont.B(:,IdStruct.IdUWRWW);
    InputMat(IdXFilt,:)         = [];

    % Meas Mat
    MeasMat      = zeros(IdStruct.NY, NX-NFilt);
    MeasMat(IdStruct.IdYGyro,:) = ModelStruct.Sens.MeasMatgyro;
    MeasMat(IdStruct.IdYAMU,:)  = ModelStruct.Sens.MeasMatAMU;
    MeasMat(IdStruct.IdYRW,:)   = ModelStruct.Sens.MeasMatRW;
    MeasMat(IdStruct.IdYSST,:)  = ModelStruct.Sens.MeasMatSST;
    MeasMat(:,RealID.IdXmod-NFilt)    = MeasMat(:,RealID.IdXmod-NFilt).*...
        (1+MeasMatCov.^0.5*randn(IdStruct.NY,2*RealID.Nmod));
    
    % Feedthrough Matrix D
    FeedThroMat = zeros(IdStruct.NY, IdStruct.NUW);
    FeedThroMat(IdStruct.IdYAMU,:) = ModelStruct.Sens.FeedThroMatAMU;

    DescMat = blkdiag(eye(RealID.NXq),M);%+NFilt
    sysMeas = dss(StateMat,InputMat,MeasMat,FeedThroMat,DescMat);
    % append sensor dynamics (anti-aliasing filter)
    sysMeas   = series(sysMeas, ModelStruct.Sens.SysSens);

    % RealStruct.RealSysDisc = c2d(sysMeas, RealStruct.SampTime, 'zoh');
    RealStruct.RealSysCont = sysMeas;
    RealStruct.RealSysContA = RealStruct.RealSysCont.A;
    RealStruct.RealSysContB = RealStruct.RealSysCont.B;
    RealStruct.RealSysContA = RealStruct.RealSysCont.C;
    RealStruct.RealSysContD = RealStruct.RealSysCont.D;

    BigA = [sysMeas.E\sysMeas.A, sysMeas.E\sysMeas.B;...
            zeros(IdStruct.NUW,RealID.NX+IdStruct.NUW)];
    AA = expm(BigA*RealStruct.SampTime);
    RealStruct.RealSysDisc = ss(AA(1:RealID.NX,1:RealID.NX),...
                                AA(1:RealID.NX,RealID.NX+1:end),...
                                sysMeas.C, ...
                                sysMeas.D, ...
                                RealStruct.SampTime);
    RealStruct.RealSysDiscA = RealStruct.RealSysDisc.A;
    RealStruct.RealSysDiscB = RealStruct.RealSysDisc.B;
    RealStruct.RealSysDiscC = RealStruct.RealSysDisc.C;
    RealStruct.RealSysDiscD = RealStruct.RealSysDisc.D;

    RealStruct.GDistJacDisc = RealStruct.RealSysDisc.B(:,IdStruct.IdUWWRW);
end

