function [ModelStruct, IdStruct, sysMeas, SysMeasDisc] = ...
                                ModelMain(ModelStruct, IdStruct, SimParam)
    file = SimParam.PostPro.plotfile;
    %% Change from FEM node coordinates to mode shape coordinates
    [ModelStruct, IdStruct] = Node2Modefunc(ModelStruct, IdStruct, SimParam);

    %% Position of vibration damping Actuator + Sensor Slots
    % A set of Nctrlpoints bundles of 3 XYZ facing RW (Reaction Wheels), and IMU 
    % sensors (XYZ Gyro + XYZ Accelerometer) are distributed along the structure. 
    % Nodes of the FEM mesh are selected.     
    [ModelStruct,IdStruct] = ActSensPlacement(ModelStruct, SimParam, IdStruct);

    %% Actuator model definition
    % Number of inputs / actuators
    % (1) X, (2) Y, (3) Z facing RWs
    [ModelStruct,IdStruct] = ActuatorModel(ModelStruct, IdStruct, SimParam);

    %% Perturbation model definition
    % Perturbations from RW jitter
    [ModelStruct,IdStruct] = PerturbationModel(ModelStruct, ...
                                                   IdStruct, SimParam);
    
    %% Sensor Models
    % There are three types of measurements in this system:
    % 1. IMU Accelerometer measurements/AMU: sense the local translational
    % acceleration
    % 2. IMU Gyroscope measurements/Gyro: sense the local angular velocity
    % Number of sensors 
    % 3. RW Encoder: measures the RW Speed
    [ModelStruct, IdStruct] = SensorModel(ModelStruct, IdStruct, SimParam);

    %% Plot Mesh
    figID = 1;
    PlotVibModes(figID, ModelStruct, 5, file)
    
    % Plot the disposition of actuators + sensors on the FEM mesh
    figID = figID + 1;
    PlotRigSysMeshActSens(figID, ModelStruct, IdStruct, file)

    %% Linearized model
    [sysMeas,IdStruct] = MergeModels(ModelStruct, IdStruct, SimParam);

    %% Get the optimal rw speed box
    figID = figID + 1;
    [ModelStruct,sysMeas] = FindOptBdwth(figID, sysMeas, ModelStruct, ...
                                                SimParam, IdStruct, file);

    %% Discretize Model
    [SysMeasDisc, IdStruct, ModelStruct] = ...
                    DiscretizeModel(sysMeas, ModelStruct, SimParam, IdStruct);

end