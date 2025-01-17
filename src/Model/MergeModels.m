function [sysMeasFilt,IdStruct] = MergeModels(ModelStruct,IdStruct, SimParam)
    % Descriptor State Space Model
    % E*xdot = A*x + B*u
    %      y = C*x + D*u
    % States: x = [q^T qdot^T]^T
    % Position states: q = eta 
    % Velocity states: qdot = [etadot^T AngVelRW^T]^T
    % Measurements: y = [GyroMeas^T AMUMeas^T AngVelRWMeas^T]^T
    % % load vars

    MeasMatAMU     = ModelStruct.Sens.MeasMatAMU;
    MeasMatRW      = ModelStruct.Sens.MeasMatRW;
    MeasMatgyro    = ModelStruct.Sens.MeasMatgyro;
    MeasMatSST     = ModelStruct.Sens.MeasMatSST;
    BInputJac      = ModelStruct.Act.BInputJac;
    GDistJac       = ModelStruct.Pert.GDistJac;
    FeedThroMatAMU = ModelStruct.Sens.FeedThroMatAMU;
    DescMat        = ModelStruct.DescMat;
    StateMat       = ModelStruct.StateMat;
    RealId         = IdStruct.RealMdl;
    % % merge models
    
    % Input Matrix B
    InputMat                             = zeros(RealId.NX, IdStruct.NUW);
    InputMat(RealId.IdXd,IdStruct.IdUWU) = BInputJac;
    InputMat(:,IdStruct.IdUWW)           = GDistJac;
    
    
    % Measurement Model
    % gyro + amu + rw speed
    % Measurement Matrix C
    %   Gyro [rad/s] | AMU [m/s^2] | AngVelRW [rad/s] | Star Tracker
    MeasMat = zeros(IdStruct.NY, RealId.NX);
    MeasMat(IdStruct.IdYGyro,:) = MeasMatgyro;
    MeasMat(IdStruct.IdYAMU,:) = MeasMatAMU;
    MeasMat(IdStruct.IdYRW,:) = MeasMatRW;
    MeasMat(IdStruct.IdYSST,:) = MeasMatSST;
    
    % Feedthrough Matrix D
    FeedThroMat = zeros(IdStruct.NY, IdStruct.NUW);
    FeedThroMat(IdStruct.IdYAMU,:) = FeedThroMatAMU;
    
    % Name of the states
    StateName = cell(1,RealId.NX);
    if RealId.Nrig
        StateName(RealId.IdXrigtq) = {'r_x', 'r_y', 'r_z'};
        StateName(RealId.IdXrigrq) = {'theta_x', 'theta_y', 'theta_z'};
        StateName(RealId.IdXrigtd) = {'v_x','v_y','v_z'};
        StateName(RealId.IdXrigrd) = {'w_x','w_y','w_z'};
    end
    StateName(RealId.IdXmodq) = cellstr(compose('eta_%d',  1:RealId.Nmod));
    StateName(RealId.IdXmodd) = cellstr(compose('etadot_%d', 1:RealId.Nmod));
    StateName(RealId.IdXRW) = cellstr(compose('AngVelRW_%d', 1:IdStruct.NU));
    
    % Name of the Measurements/Outputs
    OutputName                   = cell(1,IdStruct.NY);
    OutputName(IdStruct.IdYGyro) = cellstr(compose('gyro_%d', 1:IdStruct.NYGyro));
    OutputName(IdStruct.IdYAMU)  = cellstr(compose('accel_%d', 1:IdStruct.NYAMU));
    OutputName(IdStruct.IdYRW)   = cellstr(compose('AngVelRW_%d', 1:IdStruct.NYRW));
    OutputName(IdStruct.IdYSST)  = cellstr(compose('SST_%d', 1:IdStruct.NYSST));

    OutputName                   = cell(1,IdStruct.NY);
    OutputName(IdStruct.IdYGyro) = cellstr(compose('gyro_%d', 1:IdStruct.NYGyro));
    OutputName(IdStruct.IdYAMU)  = cellstr(compose('accel_%d', 1:IdStruct.NYAMU));
    OutputName(IdStruct.IdYRW)   = cellstr(compose('AngVelRW_%d', 1:IdStruct.NYRW));
    OutputName(IdStruct.IdYSST)  = cellstr(compose('SST_%d', 1:IdStruct.NYSST));
    
    % Name of the Actuation variables/Inputs
    InputName                 = cell(1,IdStruct.NUW);
    InputName(IdStruct.IdUWU) = cellstr(compose('U_Tz_rw_%d', 1:IdStruct.NU));
    WRWMFlag = ~isempty(IdStruct.IdUWRWMW);
    WRWFFlag = ~isempty(IdStruct.IdUWRWFW);
    NPerPerRW = 2*WRWMFlag+2*WRWMFlag;
    for iu=1:IdStruct.NU
        if WRWFFlag
            IdWFXY   = IdStruct.NU+NPerPerRW*(iu-1) + (1:2);
            InputName(IdWFXY) = compose({'W_Fx_rw_%d','W_Fy_rw_%d'}, iu);
            if WRWMFlag
                IdWMXY   = IdStruct.NU+NPerPerRW*(iu-1) + (3:4);
                InputName(IdWMXY) = compose({'W_Mx_rw_%d','W_My_rw_%d'}, iu);
            end
        else
            if WRWMFlag
                IdWMXY   = IdStruct.NU+NPerPerRW*(iu-1) + (1:2);
                InputName(IdWMXY) = compose({'W_Mx_rw_%d','W_My_rw_%d'}, iu);
            end
        end
    end
    
    % % Full State Descriptor System
    % SysFullX = dss(StateMat,InputMat,MeasMatFullX,FeedThroMatFullX,DescMat,...
    %     'StateName',StateName, 'OutputName',StateName, 'InputName',InputName);
    
    % Gyro + Accel + RW Angular Velocity Measurement Descriptor System
    sysMeas = dss(StateMat,InputMat,MeasMat,FeedThroMat,DescMat,...
        'StateName',StateName, 'OutputName',OutputName,'InputName',InputName);

    % adding sensor dynamics (anti-aliasing filter)
    sysMeasFilt = series(sysMeas, ModelStruct.Sens.SysSens);
    % system without filter
    ModelStruct.sysMeasNoFilt = sysMeas;

    %% Update Sensor indexes (after merging sensor and dynamics)
    NFiltSST     = RealId.AAFiltOrd*IdStruct.NYSST * SimParam.AASSTFlag;
    NFiltAAAMU   = RealId.AAFiltOrd*IdStruct.NYAMU * SimParam.AAAMUFlag;
    NFiltHPAMU   = RealId.HPAMUFlag*IdStruct.NYAMU;
    NFiltAMU     = NFiltAAAMU + NFiltHPAMU;
    NFiltGyr     = RealId.AAFiltOrd*IdStruct.NYGyro * SimParam.AAGyrFlag;
    NFiltIMU     = NFiltGyr + NFiltAMU;
    NAAFilt      = NFiltGyr + NFiltSST + NFiltAAAMU;
    NFilt        = NFiltGyr + NFiltSST + NFiltAMU;
    
    % update system states
    RealId.NFilt        = NFilt;
    RealId.IdXFilt      = 1:NFilt;
    RealId.NFiltHPAMU   = NFiltHPAMU;
    RealId.IdXFiltHPAMU = 1:NFiltHPAMU;
    RealId.NAAFilt      = NAAFilt;
    RealId.IdXAAFilt    = NFiltHPAMU + (1:NAAFilt);
    RealId.NFiltGyr     = NFiltGyr; 
    RealId.IdXAAFilt    = NFiltHPAMU + (1:NFiltGyr);
    RealId.NFiltAAAMU   = NFiltAAAMU; 
    RealId.IdXAAFiltAMU = NFiltHPAMU + NFiltGyr + (1:NFiltAAAMU);
    RealId.NFiltIMU     = NFiltIMU;
    RealId.IdXFiltIMU   = (1:NFiltIMU);
    RealId.NFiltSST     = NFiltSST;
    RealId.IdXFiltSST   = NFiltHPAMU + NFiltIMU + (1:NFiltSST);
    
    % total number of position states and indices within them
    RealId.IdXq     = RealId.NFilt + (1:RealId.NXq);
    RealId.IdXrigq  = RealId.NFilt + RealId.IdXqrig;
    RealId.IdXrigtq = RealId.NFilt + RealId.IdXqrigt;
    RealId.IdXrigrq = RealId.NFilt + RealId.IdXqrigr;
    RealId.IdXmodq  = RealId.NFilt + RealId.IdXqmod;
    
    % total number of velocity states and indices within them
    RealId.IdXRW    = RealId.NFilt + RealId.NXq + RealId.IdXdRW;
    RealId.IdXd     = RealId.NFilt + RealId.NXq + (1:RealId.NXd);
    RealId.IdXrigd  = RealId.NFilt + RealId.NXq + RealId.IdXdrig;
    RealId.IdXrigtd = RealId.NFilt + RealId.NXq + RealId.IdXdrigt;
    RealId.IdXrigrd = RealId.NFilt + RealId.NXq + RealId.IdXdrigr;
    RealId.IdXmodd  = RealId.NFilt + RealId.NXq + RealId.IdXdmod;
    
    % total number and Id of rigid body coordinates
    RealId.IdXmod    = [RealId.IdXmodq, RealId.IdXmodd];
    RealId.NX        = RealId.NFilt + RealId.NXq + RealId.NXd;
    IdStruct.RealMdl = RealId;
end

