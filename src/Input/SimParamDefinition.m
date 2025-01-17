function SimParam = SimParamDefinition(SimParam, ScenarioFlag)

    if strcmp(ScenarioFlag, 'Slewing')
        % Target Angle range [min, max] degrees
        SimParam.TgtSlewAngleRng = [0.5, 0.5];
        % Jitter Perturbations Flag
        SimParam.SimRWPertFlag = false; % true; % 
        % Initial state error 
        SimParam.InitAttErr = 0*((1/3600)*pi/180)^2; % [rad^2]
        SimParam.tSim     = 300; % simulation time [s]
    elseif strcmp(ScenarioFlag, 'FinePointing')
        % Target Angle range [min, max] degrees
        SimParam.TgtSlewAngleRng = 0*[0.5, 0.5];
        % Jitter Perturbations Flag
        SimParam.SimRWPertFlag = true; % false; 
        % Initial state error 
        SimParam.InitAttErr = ((1/3600)*pi/180)^2; % [rad^2]
        SimParam.tSim     = 300; % simulation time [s]
    else
        msg = ['Scenario Flag ' ScenarioFlag ' not valid.' ...
                ' Valid Options are "Slewing" and "FinePointing'];
        error(msg)
    end
    % % Simulation Configuration 
    SimParam.NMC      = 1; % number of MC runs (1 = single run)
    SimParam.OpenLoop = false; % true; % true = Open-loop simulation
    
    % True Method of propagation: ['linear', 'rk4', 'ode15s']
    SimParam.propmet  ='rk4'; % 'ode15s'; % 'linear'; %  
    SimParam.propopts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    % TODO: True Measurement Function: ['linear', 'nonlinear']
    SimParam.obsmet   = 'nonlinear'; %  'linear'; %
    % include measurement error flag
    SimParam.measerr  = true; % false; % 
    % % Slewing scenario
    % Slew Start Time
    SimParam.TSlewStart = 5;
    SimParam.raninitrw  = false;
    % SysId Model Uncertainty (unless actual sysid is done on the real model,
    % this circumvents the process and defines the resulting uncertainty for MC
    % campaigns)
    % [Ixx,Iyy,Izz,Ixy,Ixz,Iyz] % in relative value
    SimParam.InertiaCov       = 0*(1e-3)^2;
    SimParam.ModPartRotMatCov = 0*(1e-3)^2;
    SimParam.ModPartLinMatCov = 0*(1e-3)^2;
    SimParam.omegaCov         = 0*(5e-3)^2;
    SimParam.dampCov          = 0*(5e-2)^2;
    SimParam.BInputJacCov     = 0*(1e-3)^2;
    SimParam.MeasMatCov       = 0*(5e-5)^2;
end

