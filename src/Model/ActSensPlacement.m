function [ModelStruct,IdStruct] = ActSensPlacement(ModelStruct, SimParam, IdStruct)
    % % load vars
    NodeCoord = ModelStruct.NodeCoord;
    IdInner   = ModelStruct.IdInner;

    % FEM mesh nodes where controllers are placed (sensor always collocated
    % with controller)
    IdStruct.Nctrlpoints = SimParam.NCollActSens; % 3DoF RW scheme per point
    % FEM mesh nodes where additional IMUs are placed (beyond the ones for
    % controllers)
    ExtraobsPoints       = SimParam.NExtraSens;
    IdStruct.Nobspoints  = IdStruct.Nctrlpoints + ExtraobsPoints; % 3DoF IMU scheme per point

    % If user defined 
    % Select positions for the actuators and collocated sensors 
    ActNodesID  = SimParam.CollActSensNodeId;
    SensNodesID = [SimParam.CollActSensNodeId;
                    SimParam.ExtraSensNodeId];

    if ~isempty(ActNodesID)
        [~,SelActId] = find(IdInner == ActNodesID);
        SelActNodes  = NodeCoord(ActNodesID,:);
    else
        SelActId    = [];
        SelActNodes = [];
    end

    if ~isempty(SensNodesID)
        [~,SelSensId] = find(IdInner == SensNodesID);
        SelSensNodes  = NodeCoord(SensNodesID,:);
    else
        SelSensId    = [];
        SelSensNodes = [];
    end

    % Add Central rigid body actuator (if rigid body modes)
    SelActNodes  = [0,0;SelActNodes];
    SelSensNodes = [0,0;SelSensNodes];
    
    %% save vars
    % Id within interior nodes
    ModelStruct.SelActId     = SelActId;
    ModelStruct.SelSensId    = SelSensId;
    ModelStruct.SelActNodes  = SelActNodes;
    ModelStruct.SelSensNodes = SelSensNodes;
    % Id within all nodes
    ModelStruct.ActNodesID   = ActNodesID;
    ModelStruct.SensNodesID  = SensNodesID;

end
