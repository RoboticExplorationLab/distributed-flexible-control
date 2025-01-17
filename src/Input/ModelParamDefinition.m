function SimParam = ModelParamDefinition(SimParam, ModelStruct, ...
                           ActPlacementFlag, ScenarioFlag, GuidProfileFlag)
    
    %% Configuration
    % % Real plate model
    % Rayleigh damping coefficients 
    % source: www.simscale.com/knowledge-base/rayleigh-damping-coefficients/
    wn1 = 0.5; wn2 = 10; % 50; %   % lower and upper cutoff damping frequencies
    xeta1 = 0.005; xeta2 = 0.005; % damping coefficient at frequency 1 and 2
    beta = [2*pi*wn1, (1/(2*pi*wn1)); wn2*2*pi, (1/(wn2*2*pi))]\(2*[xeta1;xeta2]); 
    SimParam.beta     = beta;
    % number of modes of real model
    SimParam.Nmodes   = 25; % 15; % 5; % 10; % 50; % 395; % 100; % 200; % 
    SimParam.modeRng  = [0.1,80]; % [0.1,10]; % range of frequencies [Hz] for mode sel

    if strcmp(ActPlacementFlag,'Centralized')
        SimParam.NCollActSens = 1; % number of collocated actuators and sensors
        % number of extra sensors
        SimParam.NExtraSens   = 10; % 12; % 
        % Node of Actuators and Sensors (when user defined) - needs to be column
        % to choose nodes
        SimParam.CollActSensNodeId = [];
        SimParam.ExtraSensNodeId = [134; 68; 88; 76; 80; 84; 72;135;115;116];
        SimParam.RWperpoint    = [6]; % [3];
        SimParam.RWDirPerPoint = [eye(3);-eye(3)]; % eye(3);
        ScaleFactor = 11/2;
    elseif strcmp(ActPlacementFlag,'Distributed')
        SimParam.NCollActSens = 11; % number of collocated actuators and sensors
        % number of extra sensors
        SimParam.NExtraSens   = 0;
        SimParam.CollActSensNodeId = [134; 68; 88; 76; 80; 84; 72;135;115;116];
        SimParam.ExtraSensNodeId   = [];
        SimParam.RWperpoint    = [3,3*ones(1,10)];
        SimParam.RWDirPerPoint = [eye(3);
                                  repmat([eye(3);-eye(3)],5,1)];
        ScaleFactor = 1;
    else
        msg = ['Actuator Placement Flag ' ActPlacementFlag ' not valid.' ...
                ' Valid Options are "Centralized" and "Distributed'];
        error(msg)
    end

    % Candidate nodes for optimal placement
    SimParam.CandActNodeID     = ModelStruct.IdInner;
    % remove boundary nodes
    BoundCandIDs = any(SimParam.CandActNodeID' == ModelStruct.IdBound',1);
    SimParam.CandActNodeID(BoundCandIDs) = [];
    
    % % Model Configuration 
    % % % System Inputs
    % Dimensions of RW (X8318S) [m] 
    SimParam.RadRW       = 0.916/2; % radius of RW cylinder [m]
    SimParam.HghtRRW     = 0.41;    % height of RW cylinder [m]
    SimParam.RWInertia   = ScaleFactor*0.84; % 1.295;    % inertia of wheel [kgm^2]
    SimParam.MaxTorqueRW = ScaleFactor*0.8; % 12;       % max torque [Nm]
    % Reaction Wheel Speed/bandwidth choice
    SimParam.RWBdwthopt = 'optim'; % 'predef'; % ['predef','optim']
    % target bandwidth of optimal RW bandwidth (assuming optim)
    SimParam.TgtMinBdwt = 50/60*2*pi; % [rad/s]
    % minimum frequency (to avoid crossing in optim)
    SimParam.TgtMinFreq = 200/60*2*pi; % 15/60*2*pi; % [rad/s]
    % predefined values (centralized)
    SimParam.maxrwspeed = 28.4258;
    SimParam.minrwspeed = 23.1898;
    SimParam.nomrwspeed = 25.8078; 
    
    % % Perturbations
    % RW Jitter
    SimParam.RadMPertRW = true; % radial torque perturbation from RW
    SimParam.RadFPertRW = true; % radial force perturbation from RW
    
    % broadband noise level scalings x2 (Force and Torque)
    h1f = ScaleFactor*0.7852e-7; % [N/rpm^2]
    h1t = ScaleFactor*0.3239e-7; % [Nm/rpm^2]
    SimParam.FRWRadBBVar = (h1f*(60^2))^2; % [N/sqrt(Hz)]
    SimParam.MRWRadBBVar = (h1t*(60^2))^2; % [Nm/sqrt(Hz)]
    % First Harmonic scaling Coefficient x2 (F and M)
    rpm2rad = (2*pi)/60;
    SimParam.FRWRadFHCoeff = h1f/(rpm2rad^2); % [N/rad^2]
    SimParam.MRWRadFHCoeff = h1t/(rpm2rad^2); % [Nm/rad^2]
    % sensors
    % Gyro: 
    SimParam.NGyroperPoint = [3,2*ones(1,10)];
    SimParam.GyrDir        = [eye(3);repmat([1,0,0;0,1,0],10,1)];
    % 0.001°/√hr -> 4e-6, 40 times higher
    SimParam.GyroWNVar     = (1e-5*pi/180)^2; % 1e-3; % made up - look for reference
    SimParam.GyroBiasVar   = 0*(1e-6)^2;
    % SimParam.NAMUperPoint  = zeros(1,15); % ones(1,15);
    % SimParam.AMUDir        = repmat([0,0,1],15,1);
    SimParam.NAMUperPoint  = [3,1*ones(1,10)]; % zeros(1,15); % ones(1,15);
    SimParam.AMUDir        = [eye(3);repmat([0,0,1],10,1)];
    SimParam.NYSST         = 3; % star tracker assumed to observe orientation
    % filtering 
    SimParam.AAFiltFreq    = 80;
    SimParam.AAFiltOrd     = 2; % order of the butterworth anti-aliasing filter 
    SimParam.AAAMUFlag     = true;
    SimParam.AAGyrFlag     = false;
    SimParam.AASSTFlag     = false;
    SimParam.HPAMU         = [0,0.1]; % [true/false, cutoff freq [Hz]]
    
    % of the IMU and SST sensors
    % needs to be adjusted to go from cont to disc
    SimParam.AMUWNVar      = (2e-5)^2; 
    SimParam.AMUBiasVar    = 0*(1e-5)^2; 
    SimParam.RWWNVar       = (3e-4)^2; 
    
    % Sensor Sampling Frequency
    SimParam.SamplingFreq  = 200; % 10; % 20; % [Hz]
    SimParam.SSTSampFreq   = 10; % [Hz]
    SimParam.NStepSST      = round(SimParam.SamplingFreq / SimParam.SSTSampFreq);

    if strcmp(ScenarioFlag, 'Slewing')
        SimParam.SSTWNVar      = ((10/3600)*pi/180)^2; % [rad]
    elseif strcmp(ScenarioFlag, 'FinePointing')
        SimParam.SSTWNVar      = ((0.1/3600)*pi/180)^2; % [rad]
    else
        msg = ['Scenario Flag ' ScenarioFlag ' not valid.' ...
                ' Valid Options are "Slewing" and "FinePointing'];
        error(msg)
    end

    %% Post-processing
    SimParam.PostPro           = struct;
    SimParam.PostPro.NVib      = 5;     % Number of Vibration modes to plot
    SimParam.PostPro.att       = true;  % attitude plots
    SimParam.PostPro.SensPlots = false; % sensor/estimation plots
    SimParam.PostPro.ActPlots  = true;  % actuation plots
    SimParam.PostPro.vibmodes  = true;  % vibration modes plots
    if strcmp(ActPlacementFlag,'Centralized') % ActPlacementFlag
        if strcmp(ScenarioFlag,'Slewing') % ScenarioFlag
            if strcmp(GuidProfileFlag,'TimeOptimal') % GuidProfileFlag
                SimParam.PostPro.plotfile  = 'Results\Plots\CentSlewTO\';
                SimParam.PostPro.datafile  = 'Results\SimData\CentSlewTO.mat';
            elseif strcmp(GuidProfileFlag,'Smoothed') 
                SimParam.PostPro.plotfile  = 'Results\Plots\CentSlewSmoothed\';
                SimParam.PostPro.datafile  = 'Results\SimData\CentSlewSmoothed.mat';
            else
                msg = ['Guidance Profile Flag ' GuidProfileFlag ' not valid.' ...
                ' Valid Options are "TimeOptimal" and "Smoothed'];
                error(msg)
            end
        elseif strcmp(ScenarioFlag,'FinePointing') % ScenarioFlag
            SimParam.PostPro.plotfile  = 'Results\Plots\CentFinePointing\';
            SimParam.PostPro.datafile  = 'Results\SimData\CentFinePointing.mat';
        else
            msg = ['Scenario Flag ' ScenarioFlag ' not valid.' ...
                ' Valid Options are "Slewing" and "FinePointing'];
            error(msg)
        end
    elseif strcmp(ActPlacementFlag,'Distributed')
        if strcmp(ScenarioFlag,'Slewing') % ScenarioFlag
            if strcmp(GuidProfileFlag,'TimeOptimal') % GuidProfileFlag
                SimParam.PostPro.plotfile  = 'Results\Plots\DistSlewTO\';
                SimParam.PostPro.datafile  = 'Results\SimData\DistSlewTO.mat';
            elseif strcmp(GuidProfileFlag,'Smoothed') 
                SimParam.PostPro.plotfile  = 'Results\Plots\DistSlewSmoothed\';
                SimParam.PostPro.datafile  = 'Results\SimData\DistSlewSmoothed.mat';
            else
                msg = ['Guidance Profile Flag ' GuidProfileFlag ' not valid.' ...
                ' Valid Options are "TimeOptimal" and "Smoothed'];
                error(msg)
            end
        elseif strcmp(ScenarioFlag,'FinePointing') % ScenarioFlag
            SimParam.PostPro.plotfile  = 'Results\Plots\DistFinePointing\';
            SimParam.PostPro.datafile  = 'Results\SimData\DistFinePointing.mat';
        else
            msg = ['Scenario Flag ' ScenarioFlag ' not valid.' ...
                ' Valid Options are "Slewing" and "FinePointing'];
            error(msg)
        end
    else
        msg = ['Actuator Placement Flag ' ActPlacementFlag ' not valid.' ...
                ' Valid Options are "Centralized" and "Distributed'];
        error(msg)
    end

end

