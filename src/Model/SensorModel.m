function [ModelStruct, IdStruct] = SensorModel(ModelStruct, IdStruct, SimParam)
    % % load vars
    % SelSensNodes SelSensId Node2Mode
    SelSensNodes = ModelStruct.SelSensNodes;
    SelSensId    = ModelStruct.SelSensId;
    Node2Mode    = ModelStruct.Node2Mode;
    ConnectMat   = ModelStruct.ConnectMat;
    NodeCoord    = ModelStruct.NodeCoord;
    StateMat     = ModelStruct.StateMat;
    BInputJac    = ModelStruct.Act.BInputJac;
    DescMat      = ModelStruct.DescMat;
    StateMatEA   = DescMat\StateMat;
    GDistJac     = ModelStruct.Pert.GDistJac;
    RealId       = IdStruct.RealMdl;
    IdXd         = RealId.IdXd;
    BInputJacEB  = DescMat(IdXd,IdXd)\[BInputJac,GDistJac(IdXd,:)];

    Nfuncs   = ModelStruct.Nfuncs;
    IdBound    = ModelStruct.IdBound*3 + (-2:0);
    IdBound    = IdBound(:);

    % % Sensor Indexes
    IdStruct.NGyroperPoint = SimParam.NGyroperPoint;
    IdStruct.NAMUperPoint  = SimParam.NAMUperPoint;
    IdStruct.NYSST         = SimParam.NYSST;
    IdStruct.NYGyro  = sum(IdStruct.NGyroperPoint);
    IdStruct.NYAMU   = sum(IdStruct.NAMUperPoint);
    IdStruct.NYIMU   = IdStruct.NYGyro + IdStruct.NYAMU;
    IdStruct.NYRW    = RealId.NXRW;
    IdStruct.NY      = IdStruct.NYIMU + IdStruct.NYRW + IdStruct.NYSST;
    
    IdStruct.IdYGyro = 1:IdStruct.NYGyro;
    IdStruct.IdYAMU  = IdStruct.NYGyro + (1:IdStruct.NYAMU);
    IdStruct.IdYIMU  = [IdStruct.IdYGyro, IdStruct.IdYAMU];
    IdStruct.IdYRW   = IdStruct.NYIMU + (1:IdStruct.NYRW);
    IdStruct.IdYSST  = IdStruct.NYIMU + IdStruct.NYRW + (1:IdStruct.NYSST);

    % Position of the IMUs
    IMUPos          = [SelSensNodes, zeros(IdStruct.Nobspoints,1)];
    % Orientation of the axis of the IMU
    GyrDir          = SimParam.GyrDir;
    AMUDir          = SimParam.AMUDir;

    % Input Jacobian link between gyro and nodal rotation around
    % rotation axis
    Lgyr = zeros(IdStruct.NYGyro, RealId.Nmod);
    % Input Jacobian link between gyro and nodal rotation
    Ngyr = zeros(IdStruct.NYGyro,3, RealId.Nmod);
    % Input Jacobian link between accel and nodal displacement along
    % displacement axis
    Lamu = zeros(IdStruct.NYAMU, RealId.Nmod);
    % Input Jacobian link between accel and nodal displacement
    Namu = zeros(IdStruct.NYAMU, RealId.Nmod,3);

    % Measurement Matrix of the Gyro
    MeasMatgyro     = zeros(IdStruct.NYGyro, RealId.NX);
    MeasMatAMU      = zeros(IdStruct.NYAMU, RealId.NX);
    MeasMatAMUaccel = zeros(IdStruct.NYAMU, RealId.NX);
    FeedThroMatAMU  = zeros(IdStruct.NYAMU,IdStruct.NUW);
    MeasMatSST      = zeros(IdStruct.NYSST, RealId.NX);
    MeasMatSST(:, RealId.IdXrigrq) = 4*eye(3);

    for i=1:IdStruct.Nobspoints
        % % Gyro meas = dPhidx * etadot
        % % Gyro meas = dPhidx * etadot
        if IdStruct.NGyroperPoint
            IdGyrI = sum(IdStruct.NGyroperPoint(1:i-1))+1:sum(IdStruct.NGyroperPoint(1:i));
            MeasMatgyro(IdGyrI,RealId.IdXrigrd) = GyrDir(IdGyrI,:);
            if i>1
                % femN      = FemN(IMUPos(i,:),Nfuncs,ConnectMat,NodeCoord, ...
                %                                         IdBound)*Node2Mode;
                % Crot      = [femN(2,:);femN(3,:);zeros(1,IdStruct.LinMdl.Nmod)];
                % MeasMatgyro(IdGyrI,IdStruct.LinMdl.IdXmodd) = GyrDir(IdGyrI,:)*Crot;
                IdGyr = SelSensId(i-1)*3 + (-1:0); % (-2:-1);
                Ngyr(i,1:2,:) = Node2Mode(IdGyr,:);
                Lgyr(IdGyrI,:) = tensorprod(GyrDir(IdGyrI,:),Ngyr(i,:,:),2);
                MeasMatgyro(IdGyrI, RealId.IdXmodd) = Lgyr(IdGyrI,:);
            end
        end

        % AMU meas  = - r x omegadot + Phi * etaddot [m/s^2]
        if IdStruct.NAMUperPoint
            IdAMUI = sum(IdStruct.NAMUperPoint(1:i-1))+1:sum(IdStruct.NAMUperPoint(1:i));
            % femN      = FemN(IMUPos(i,:),Nfuncs,ConnectMat,NodeCoord, ...
            %                                         IdBound)*Node2Mode;
            % femNr = [zeros(2, RealId.Nmod);femN(1,:)];

            if i>1
                IdAct = SelSensId(i-1)*3 + -2; % (-2:-1);
                femNr = [zeros(2,IdStruct.RealMdl.Nmod); Node2Mode(IdAct,:)];
            else
                femNr = zeros(3,IdStruct.RealMdl.Nmod);
            end

            MeasMatAMUaccel(IdAMUI,:) = ...
                       AMUDir(IdAMUI,:)*[zeros(3, RealId.NXq), ...
                                                       eye(3), ...
                                           -skew(IMUPos(i,:)), ...
                                                        femNr, ...
                                        zeros(3, RealId.NXRW)];

            MeasMatAMU(IdAMUI,:)     = MeasMatAMUaccel(IdAMUI,:)*StateMatEA;
            FeedThroMatAMU(IdAMUI,:) = ...
            MeasMatAMUaccel(IdAMUI,RealId.IdXd)*BInputJacEB;
        end
    end
    % inertial approximation
    MeasMatRW                 = zeros(IdStruct.NYRW,RealId.NX);
    MeasMatRW(:,RealId.IdXRW) = eye(IdStruct.NYRW);
    % more accurate 
    MeasMatRW(:,RealId.IdXrigrd) = -ModelStruct.Act.RWDir;
    MeasMatRW(:,RealId.IdXmodd)  = -ModelStruct.Act.Lrw';
    MeasMatRW                    = MeasMatRW/(2*pi); % changing to rev/s

    %% Anti-aliasing filter on analog measurement signal
    % sampling frequency
    SamplingFreq   = SimParam.SamplingFreq; % [Hz]
    % sampling time
    SamplingTime   = 1/SamplingFreq; % [s]

    % anti-aliasing analog signal (assuming sensors have anti-aliasing)
    % without anti-aliasing ERA may show weird results
    % lowpass = tf(1,[(SamplingTime/pi)^2, 2*SamplingTime/(sqrt(2)*pi), 1]);
    % order of lowpass anti-aliasing filter
    RealId.AAFiltOrd = SimParam.AAFiltOrd;
    AAFiltOrd        = SimParam.AAFiltOrd;
    % requires 2024a
    % [b,a] = butter(AAFiltOrd,SamplingFreq*pi,'s');
    % alternative if older matlab version and order 2
    if RealId.AAFiltOrd == 0
        % SysMeasRedlp = SysMeasRed;
        SysSens = ss([],[],[],eye(IdStruct.NY));
    else
        AAAMUFlag = SimParam.AAAMUFlag;
        AAGyrFlag = SimParam.AAGyrFlag;
        AASSTFlag = SimParam.AASSTFlag;
        AAFreq    = SimParam.AAFiltFreq;
        if RealId.AAFiltOrd == 2
            b = [0,0,(AAFreq*2*pi)^2];
            a = [1,sqrt(2)*AAFreq*2*pi,(AAFreq*2*pi)^2];
        else
            [b,a] = butter(AAFiltOrd,AAFreq*2*pi,'s');
        end
        [a,b,c,d] = tf2ss(b,a);
        filterstatenames = [];
        if AAGyrFlag
            lpssgyr = ss(ss(a,b,c,d));
            gyroids = cat(2,repelem(IdStruct.IdYGyro,AAFiltOrd)',...
                        repmat(1:AAFiltOrd,1,IdStruct.NYGyro)');
            filterstatenames = [filterstatenames; ...
                cellstr(compose('gyro_%d_aalp%d', gyroids(:,1),gyroids(:,2)))];
        else
            lpssgyr = 1;
        end
        if AAAMUFlag
            lpssamu = ss(ss(a,b,c,d));
            amuids = cat(2,repelem(IdStruct.IdYAMU,AAFiltOrd)',...
                       repmat(1:AAFiltOrd,1,IdStruct.NYAMU)');
            filterstatenames = [filterstatenames; ...
                cellstr(compose('amu_%d_aalp%d', amuids(:,1),amuids(:,2)))];
        else
            lpssamu = 1;
        end
        if AASSTFlag
            lpsssst = ss(ss(a,b,c,d));
            sstids = cat(2,repelem(IdStruct.IdYSST,AAFiltOrd)',...
                   repmat(1:AAFiltOrd,1,IdStruct.NYSST)');
            filterstatenames = [filterstatenames; ...
                cellstr(compose('sst_%d_aalp%d', sstids(:,1),sstids(:,2)))];
        else
            lpsssst = 1;
        end
        
        % [a,b,c,d]=tf2ss(1,[(SamplingTime/(pi))^2, 2*SamplingTime/(sqrt(2)*pi), 1]);
        SysSens = blkdiag(lpssgyr*eye(IdStruct.NYGyro),...
                          lpssamu*eye(IdStruct.NYAMU),...
                                     eye(IdStruct.NU), ...
                          lpsssst*eye(IdStruct.NYSST));
        SysSens.StateName = filterstatenames;
    end

    %% Additional Filtering
    % AMU HP filtering (removing rigid body dyn)
    RealId.HPAMUFlag = SimParam.HPAMU(1);
    if SimParam.HPAMU(1) == 1
        b = [0,1];
        a = [1,SimParam.HPAMU(2)*2*pi];
        [a,b,c,d] = tf2ss(b,a);
        hpss = ss(ss(a,b,c,d));
        % append filter to syssens
        SysHPAMU = blkdiag(eye(IdStruct.NYGyro),...
                           hpss*eye(IdStruct.NYAMU),...
                           eye(IdStruct.NU+IdStruct.NYSST));
        % SysHPAMU = hpss*eye(IdStruct.NYAMU);
        SysHPAMU.StateName = cellstr(compose('amu_%d_hp', 1:IdStruct.NYAMU));
        SysSens = series(SysSens,SysHPAMU);
        % SysSens = series(SysSens,SysHPAMU,IdStruct.IdYAMU,1:IdStruct.NYAMU);
    end
    
    %% Saving results
    ModelStruct.Sens.SamplingTime   = SamplingTime;
    ModelStruct.Sens.SysSens        = SysSens;
    ModelStruct.Sens.SysSensA       = SysSens.A;
    ModelStruct.Sens.SysSensB       = SysSens.B;
    ModelStruct.Sens.SysSensC       = SysSens.C;
    ModelStruct.Sens.SysSensD       = SysSens.D;
    % Initialize measurement errors
    % ModelStruct.Sens
    ModelStruct.Sens.GyroBias       = SimParam.GyroBiasVar*randn(IdStruct.NYGyro,1);
    ModelStruct.Sens.AMUBias        = SimParam.AMUBiasVar*randn(IdStruct.NYAMU,1);
    ModelStruct.Sens.GyroWNVar      = SimParam.GyroWNVar;
    ModelStruct.Sens.AMUWNVar       = SimParam.AMUWNVar;
    ModelStruct.Sens.SSTWNVar       = SimParam.SSTWNVar;
    ModelStruct.Sens.RWWNVar        = SimParam.RWWNVar;

    % % Save vars
    ModelStruct.Sens.MeasMatgyro     = MeasMatgyro;
    ModelStruct.Sens.MeasMatAMU      = MeasMatAMU;
    ModelStruct.Sens.MeasMatRW       = MeasMatRW;
    ModelStruct.Sens.MeasMatAMUaccel = MeasMatAMUaccel;
    ModelStruct.Sens.FeedThroMatAMU  = FeedThroMatAMU;
    ModelStruct.Sens.MeasMatSST      = MeasMatSST;

    % Position Vars
    ModelStruct.Sens.AMUDir          = AMUDir;
    ModelStruct.Sens.GyrDir          = GyrDir;
    ModelStruct.Sens.AMUPos          = [repelem(IMUPos,IdStruct.NAMUperPoint,1)];
    ModelStruct.Sens.GyrPos          = [repelem(IMUPos,IdStruct.NGyroperPoint,1)];

    % Nodal Jacobian Vars
    ModelStruct.Sens.Lgyr            = Lgyr;
    ModelStruct.Sens.Ngyr            = Ngyr;
    IdStruct.RealMdl                 = RealId;
end