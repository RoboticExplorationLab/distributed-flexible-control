clear all;
close all;
clc;

%% Input Structure Configuration definition
SimParam0 = struct;
SimParam0.fem_file_name = 'InputData/AntennaMCKScen4.mat';

%% load FEM files
TableFile = load(SimParam0.fem_file_name);

%% Input pre-processing
[ModelStruct0,IdStruct0] = InputPrePro(TableFile, SimParam0);

%% Scenario Flags
ActPlacementFlagList = {'Centralized', 'Distributed'}; % 
ScenarioFlag     = 'FinePointing'; % 'Slewing'; % 
GuidProfileFlag  = 'TimeOptimal'; % 'Smoothed'; % only relevant for Slewing

for i=1:2
    %% Configuration
    ActPlacementFlag = ActPlacementFlagList{i};
    SimParam = ModelParamDefinition(SimParam0, ModelStruct0, ActPlacementFlag, ...
                                                ScenarioFlag, GuidProfileFlag);
    
    %% Real-World Model Definition
    close all;
    [ModelStruct, IdStruct, RealSysCont, RealSysDisc] = ModelMain(ModelStruct0,...
                                                    IdStruct0, SimParam);
    
    %% Observability analysis
    [TKD, IdKD] = ObsvCtrbAnalysis(RealSysDisc, ModelStruct, IdStruct, SimParam);
    
    %% Controller Design
    close all;
    SimParam = OBSWParamDefinition(SimParam, RealSysDisc, ModelStruct, ...
                                   IdStruct, TKD, IdKD, ActPlacementFlag, ...
                                   ScenarioFlag, GuidProfileFlag);
    
    [OBSWStruct, ModelStruct, IdStruct] = ControllerDesign(RealSysDisc, ...
                                    IdStruct, ModelStruct, SimParam);
    
    %% Simulation Configuration
    close all;
    SimParam = SimParamDefinition(SimParam, ScenarioFlag);
    
    %% Simulation
    SimResMC = struct;
    SimResMC.SingleSimRes = cell(1,SimParam.NMC);
    for iMC = 1:SimParam.NMC
        % Parameter dispersion (bias, model error) [No dispersion]
        RealStruct = MCDispersion(ModelStruct, RealSysCont, IdStruct, ...
                                                                    SimParam);
    
        % Run Sim
        SimResMC.SingleSimRes{iMC} = RunSingleSim(ModelStruct,OBSWStruct, ...
                                              RealStruct, IdStruct, SimParam);
    end
    
    %% Simulation Post-Processing
    PostProMain(SimResMC.SingleSimRes{1}, SimParam, IdStruct, ModelStruct)
    
    %% save simulation
    save(SimParam.PostPro.datafile,'SimResMC','RealSysDisc','SimParam',...
                                        'IdStruct','ModelStruct','OBSWStruct')

    %% Clear vars
    clear SimResMC RealSysDisc SimParam IdStruct ModelStruct OBSWStruct
end

%% Comparison Post-Processing
PostProComparison(ActPlacementFlagList,ScenarioFlag,GuidProfileFlag);
