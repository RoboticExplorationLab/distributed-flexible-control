function PlotVibModes(figID, ModelStruct, Nmodes, file)
    K_b        = ModelStruct.StiffMatModes;
    Phi        = ModelStruct.Node2Mode;
    ConnectMat = ModelStruct.ConnectMat;
    NodeCoord  = ModelStruct.NodeCoord;
    IdSelect   = ModelStruct.IdBound;
    
    NDOFs = size(K_b,1);
    [MinFreqs,IdFreq] = mink(diag(K_b)',NDOFs,2);
    Phi = Phi(:,IdFreq);
    MinFreqs = (MinFreqs.^0.5)./(2*pi);
    
    % Plot the initial surface
    x = NodeCoord(:,1);
    y = NodeCoord(:,2);
    z = zeros(size(x));
    map = [  0,   0, 0.5;
           0.5, 0.5, 0.5;
           0.5,   0,  0];
    n = size(map,1);
    Nint = 20;
    map_interp = interp1(1:n,map,linspace(1,n,Nint*(n-1)+1));
    f1 = figure(figID);
    tiledlayout(1,Nmodes)
    for j= 1:Nmodes
        % subplot(1,Nmodes,j)
        nexttile
        w = squeeze(Phi(1:3:end,j));
        z(setdiff(1:end,IdSelect)) = w; % sign(w(1))*w;
        % scale z
        z = z / max(abs(z));
        % z = Phif(1:3:end,end);
        trisurf(ConnectMat,x,y,z,'EdgeColor','interp','FaceColor','interp',...
                                                            'FaceAlpha',0.5);
        
        colormap(map_interp) % 'gray')
        % Set axis labels and title
        xlabel('X [m]','fontsize',10);
        if j == 1
        ylabel('Y [m]','fontsize',10);
        end
        zlabel('Z-axis');
        title(['$$\eta_' num2str(j) '$$'],'Interpreter','latex','fontsize',14);
        % title(['$$\eta_' num2str(j) ', \omega_' num2str(j) '$$=' num2str(MinFreqs(j)) ' Hz'],...
        %                                     'Interpreter','latex','fontsize',14);
        view(0,90)
    
        % Iterate over time (adjust the number of iterations as needed)
        axis([min(NodeCoord(:,1))-0.5, max(NodeCoord(:,1))+0.5, ...
              min(NodeCoord(:,2))-0.5, max(NodeCoord(:,2))+0.5, ...
                      min(z,[],'all')-0.1*abs(min(z,[],'all')), ...
                      max(z,[],'all')+0.1*abs(max(z,[],'all'))])
        axis manual
        clim([-1 1])
        hold on;
        % Plot the rigid Platform    
        theta = linspace(pi/2, 15*pi/6, 7);
        fill3(2*cos(theta), 2*sin(theta), 0.1*ones(1,7), 'k', 'FaceAlpha',0.5);
        grid off
        box on
        ax = gca;
        ax.FontSize = 10; 
    end
    
    c = colorbar;
    c.Label.String = 'Z Displacement';
    c.Label.FontSize = 10;
    savefig(f1,[file, 'VibModesSC.fig']);
    saveas(f1,[file, 'VibModesSC.png']);
    addpath 'C:\Users\prca2\Documents\MATLAB_Code\matlab2tikz\src'
    matlab2tikz([file, 'VibModesSC.tex'], 'standalone', true);

end

