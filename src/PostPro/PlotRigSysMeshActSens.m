function PlotRigSysMeshActSens(figID, ModelStruct, IdStruct, file)
    %PLOTRIGSYSMESHACTSENS Summary of this function goes here
    %   Detailed explanation goes here
    ConnectMat   = ModelStruct.ConnectMat;
    NodeCoord    = ModelStruct.NodeCoord;
    SelActNodes  = ModelStruct.SelActNodes;
    SelSensNodes = ModelStruct.SelSensNodes;
    NU           = IdStruct.NU;
    Nobspoints   = IdStruct.Nobspoints;
    cp = parula(3);
    dp = autumn(3);
    % [TODO]: 
    f1 = figure(figID);
    ax=gca;
    ax.FontSize = 10;
    h = cell(4+3*NU+Nobspoints,1);
    % Plot the mesh
    
    h{1} = trimesh(ConnectMat,NodeCoord(:,1),NodeCoord(:,2),...
                zeros(size(NodeCoord(:,2))),zeros(size(NodeCoord(:,2))),...
                'EdgeColor',[0.5,0.5,0.5],'FaceColor',[0.5,0.5,0.5], ...
                'FaceAlpha',0.2,'LineWidth',1.25);
    colormap([0.5 0.5 0.5; cp(2,:); dp(2,:)]);
    hold on;
    % Plot the rigid Platform    
    theta = linspace(pi/2, 15*pi/6, 7);
    x = 2*cos(theta);
    y = 2*sin(theta);
    X = [x;x];
    Y = [y;y];
    Z = [-1*ones(1,7); 2*ones(1,7)];
    h{2} = mesh(X,Y,Z,'facecolor','k', 'FaceAlpha',0.5);
    h{3} = fill3(X(1,:), Y(1,:), Z(1,:), 'k', 'FaceAlpha',0.5);
    h{4} = fill3(X(2,:), Y(2,:), Z(2,:), 'k', 'FaceAlpha',0.5);

    % Plot the Actuators
    for ia = 1:NU
        [Xrw0,Yrw0,Zrw0] = cylinder(ModelStruct.Act.RadRW(ia));
        Zrw0 = Zrw0*ModelStruct.Act.HghtRRW(1)-(ModelStruct.Act.HghtRRW(ia)/2);
        % Compute the rotation matrix
        a = [0, 0, 1];
        b = ModelStruct.Act.RWDir(ia,:);
        v = cross(a, b);
        c = dot(a, b);
        v_cross = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
        v_cross_square = v_cross * v_cross;
        rotationMatrix = eye(3) + v_cross + v_cross_square ./ (1 + c);
       
        % get points at the two rings and rotate them separately:
        positionOld1 = [Xrw0(1,:)',Yrw0(1,:)',Zrw0(1,:)'];
        positionOld2 = [Xrw0(2,:)',Yrw0(2,:)',Zrw0(2,:)'];
        positionNew1 = positionOld1*rotationMatrix;
        positionNew2 = positionOld2*rotationMatrix;
        
        % reassemble the two sets of points into X Y Z format:
        Xrw = [positionNew1(:,1),positionNew2(:,1)];
        Yrw = [positionNew1(:,2),positionNew2(:,2)];
        Zrw = [positionNew1(:,3),positionNew2(:,3)];

        % translate 
        dpos = b*(ModelStruct.Act.RadRW(ia)+(ModelStruct.Act.HghtRRW(ia)/2));
        Xrw = Xrw + ModelStruct.Act.RWPos(ia,1) + dpos(1);
        Yrw = Yrw + ModelStruct.Act.RWPos(ia,2) + dpos(2);
        Zrw = Zrw + dpos(3);
        % % rotate
        h{4+3*(ia-1)+1} = surf(Xrw,Yrw,Zrw,2*ones(size(Zrw)), ...
                                    'EdgeColor','none', 'FaceAlpha',0.8);
        h{4+3*(ia-1)+2} = fill3(Xrw(:,1), Yrw(:,1), Zrw(:,1), ...
                                2*ones(size(Zrw(:,1))), 'FaceAlpha',0.8);
        h{4+3*(ia-1)+3} = fill3(Xrw(:,2), Yrw(:,2), Zrw(:,2), ...
                                2*ones(size(Zrw(:,2))), 'FaceAlpha',0.8);
    end
    hold on;
    % Plot the IMUs
    limu = 0.20;
    coord = [-limu, -limu,       0;
              limu, -limu,       0;
              limu,  limu,       0;
             -limu,  limu,       0;
             -limu, -limu,  2*limu;
              limu, -limu,  2*limu;
              limu,  limu,  2*limu;
             -limu,  limu,  2*limu;];
    idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';

    for io = 1:Nobspoints
        xc = SelSensNodes(io,1)+coord(:,1);
        yc = SelSensNodes(io,2)+coord(:,2);
        zc = coord(:,3);
        h{4+3*NU+io} = patch('XData', xc(idx), 'YData', yc(idx), 'ZData', zc(idx), ...
            'FaceColor',cp(2,:), 'EdgeColor','k');
    end
    hold on;

    % Plot the SSTs (?)
    coord = [-limu, -limu,       0;
              limu, -limu,       0;
              limu,  limu,       0;
             -limu,  limu,       0;
             -limu, -limu,  2*limu;
              limu, -limu,  2*limu;
              limu,  limu,  2*limu;
             -limu,  limu,  2*limu;];
    idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';
    xc = coord(:,1);
    yc = coord(:,2);
    zc = coord(:,3)+2;
    h{5+3*NU+Nobspoints} = patch('XData', xc(idx), 'YData', yc(idx), ...
                            'ZData', zc(idx),'facecolor',[0.5,0,0]); 
    % , 'facealpha', 0.1);

    % title('FEM Mesh + HW placement')
    legend([h{1}(1) h{3}(1) h{6} h{5+3*NU} h{5+3*NU+Nobspoints}],...
        {'FEM Mesh','Platform','RW','IMU','SST'}, 'Location','northeast', 'FontSize',10)
    axis equal
    xlim([min(NodeCoord(:,1))-1, max(NodeCoord(:,1))+1])
    ylim([min(NodeCoord(:,2))-1, max(NodeCoord(:,2))+1])
    zlim([-1, 2+1])
    xlabel('X [m]', 'FontSize',10)
    ylabel('Y [m]', 'FontSize',10)
    zlabel('Z [m]', 'FontSize',10)
    % grid on;
    box off
    grid off;
    view([45 30])
    savefig(f1,[file, 'MeshActSens.fig']);
    saveas(f1,[file, 'MeshActSens.png']);
end