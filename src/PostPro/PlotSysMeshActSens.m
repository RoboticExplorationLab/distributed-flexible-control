function PlotSysMeshActSens(figID,ModelStruct, file)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    % % load var
    ConnectMat = ModelStruct.ConnectMat;
    NodeCoord = ModelStruct.NodeCoord;
    SelActNodes = ModelStruct.SelActNodes;
    SelSensNodes = ModelStruct.SelSensNodes;
    c = parula(3);
    d = autumn(3);
    % % plot Plot
    f1 = figure(figID);
    h = cell(6,1);
    h{1} = quadmesh(ConnectMat,NodeCoord(:,1),NodeCoord(:,2),NodeCoord(:,3),...
                    zeros(size(NodeCoord(:,3))),'FaceColor','interp',...
                                       'FaceAlpha',0.2,'LineWidth',1.25); 
    colormap([0.5 0.5 0.5; c(2,:); d(2,:)]);
    % ,'Color','b','LineWidth',1.25);
    hold on;
    [Xrw0,Yrw0,Zrw0] = cylinder(ModelStruct.Act.RadRW);
    Zrw0 = Zrw0*ModelStruct.Act.HghtRRW;
    %translate cylinder to new location
    Xrw = Xrw0 + SelActNodes(:,1); 
    Yrw = Yrw0 + SelActNodes(:,2); 
    Zrw = Zrw0 + SelActNodes(:,3); 

    h{2} = surf(Xrw,Yrw,Zrw,0.5*ones(size(Zrw)), ...
                                'EdgeColor','none', 'FaceAlpha',0.8);
    h{3} = fill3(Xrw(1,:), Yrw(1,:), Zrw(1,:), ...
                                0.5*ones(size(Zrw(1,:))), 'FaceAlpha',0.8);
    h{4} = fill3(Xrw(2,:), Yrw(2,:), Zrw(2,:), ...
                                0.5*ones(size(Zrw(2,:))), 'FaceAlpha',0.8);


    % h{2} = plot3(SelActNodes(:,1),SelActNodes(:,2),SelActNodes(:,3),'ko',...
    %                                          'MarkerSize',12,'LineWidth',1.5);
    hold on;
    
    limu = 0.010;
    coord = [-limu, -limu, -2*limu;
              limu, -limu, -2*limu;
              limu,  limu, -2*limu;
             -limu,  limu, -2*limu;
             -limu, -limu,       0;
              limu, -limu,       0;
              limu,  limu,       0;
             -limu,  limu,       0;];
    idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';
    xc = SelSensNodes(:,1)+coord(:,1);
    yc = SelSensNodes(:,2)+coord(:,2);
    zc = SelSensNodes(:,3)+coord(:,3);
    h{3} = patch('XData', xc(idx), 'YData', yc(idx), 'ZData', zc(idx), ...
        'FaceColor',c(2,:), 'EdgeColor','k'); % , 'facealpha', 0.1);
    % h{3} = plot3(SelSensNodes(:,1),SelSensNodes(:,2),SelSensNodes(:,3),...
    %                              'rx', 'MarkerSize',12,'LineWidth',1.5);
    hold on;
    coord = [-0.1, -0.1, -0.1;
                0, -0.1, -0.1;
                0,  0.1, -0.1;
             -0.1,  0.1, -0.1;
             -0.1, -0.1,  0.5;
                0, -0.1,  0.5;
                0,  0.1,  0.5;
             -0.1,  0.1,  0.5;];
    idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';
    xc = coord(:,1);
    yc = coord(:,2);
    zc = coord(:,3);
    h{4} = patch(xc(idx), yc(idx), zc(idx), 'k'); % , 'facealpha', 0.1);
    h{5} = line(NaN,NaN,'LineWidth',4,'LineStyle','-','Color',d(2,:));
    h{6} = line(NaN,NaN,'LineWidth',4,'LineStyle','-','Color',c(2,:));
    title('FEM Mesh + HW Placement')
    legend([h{1}(1) h{5} h{6}],{'FEM Mesh','RW','IMU'},'Location','northeast')
    axis equal
    xlim([min(NodeCoord(:,1))-0.1, max(NodeCoord(:,1))+0.1])
    ylim([min(NodeCoord(:,2))-0.1, max(NodeCoord(:,2))+0.1])
    zlim([min(NodeCoord(:,3))-0.1, max(NodeCoord(:,3))+0.1])
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    grid
    view([45 10])
    % TODO: 
    % 1. change shape of RW actuator to make it clear which way it will be
    %       spinning
    % 2.  

    savefig(f1,[file,'MeshActSens2DPlate.fig']);
    saveas(f1,[file, 'MeshActSens2DPlate.png']);
    % close(f1)
end

