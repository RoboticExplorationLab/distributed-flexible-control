function PostProComparison(ActPlacementFlagList,ScenarioFlag,GuidProfileFlag)
    close all;
    clc;

    %% Load results
    N = length(ActPlacementFlagList);
    InputFiles = cell(1,N);
    OutputFile = 'Results\Plots\';
    ScenNames = cell(1,N);
    for i=1:N
        ActPlacementFlag = ActPlacementFlagList(i);
        % slewing vs fine pointing
        if strcmp(ScenarioFlag,'Slewing')
            if strcmp(GuidProfileFlag,'TimeOptimal')
                if i==1
                    OutputFile = [OutputFile 'SlewTO\'];
                end
                if strcmp(ActPlacementFlag,'Centralized')
                    ScenNames{i}  = 'Cent';
                    InputFiles{i} = 'Results\SimData\CentSlewTO.mat';
                elseif strcmp(ActPlacementFlag,'Distributed')
                    ScenNames{i} = 'Dist';
                    InputFiles{i} = 'Results\SimData\DistSlewTO.mat';
                else
                    msg = ['Actuator Placement Flag ' ActPlacementFlag ' not valid.' ...
                            ' Valid Options are "Centralized" and "Distributed'];
                    error(msg)
                end
            elseif strcmp(GuidProfileFlag,'Smoothed') 
                if i==1
                    OutputFile = [OutputFile 'SlewSmoothed\'];
                end
                if strcmp(ActPlacementFlag,'Centralized')
                    ScenNames{i}  = 'Cent';
                    InputFiles{i} = 'Results\SimData\CentSlewSmoothed.mat';
                elseif strcmp(ActPlacementFlag,'Distributed')
                    ScenNames{i}  = 'Dist';
                    InputFiles{i} = 'Results\SimData\DistSlewSmoothed.mat';
                else
                    msg = ['Actuator Placement Flag ' ActPlacementFlag ' not valid.' ...
                            ' Valid Options are "Centralized" and "Distributed'];
                    error(msg)
                end
            else
                msg = ['Guidance Profile Flag ' GuidProfileFlag ' not valid.' ...
                ' Valid Options are "TimeOptimal" and "Smoothed'];
                error(msg)
            end
        elseif strcmp(ScenarioFlag,'FinePointing') 
            if i==1
                OutputFile = [OutputFile 'FinePointing\'];
            end
            if strcmp(ActPlacementFlag,'Centralized')
                ScenNames{i}  = 'Cent';
                InputFiles{i} = 'Results\SimData\CentFinePointing.mat';
            elseif strcmp(ActPlacementFlag,'Distributed')
                ScenNames{i}  = 'Dist';
                InputFiles{i} = 'Results\SimData\DistFinePointing.mat';
            else
                msg = ['Actuator Placement Flag ' ActPlacementFlag ' not valid.' ...
                        ' Valid Options are "Centralized" and "Distributed'];
                error(msg)
            end
        else
            msg = ['Scenario Flag ' ScenarioFlag ' not valid.' ...
                ' Valid Options are "Slewing" and "FinePointing'];
            error(msg)
        end
    end
    
    resdata    = cell(N,1);
    for i=1:N
        resdata{i} = load(InputFiles{i});
        for fn = fieldnames(resdata{i}.SimResMC.SingleSimRes{1})'
           resdata{i}.(fn{1}) = resdata{i}.SimResMC.SingleSimRes{1}.(fn{1});
        end
        for fn = fieldnames(resdata{i}.ModelStruct)'
           resdata{i}.(fn{1}) = resdata{i}.ModelStruct.(fn{1});
        end
        for fn = fieldnames(resdata{i}.ModelStruct.Act)'
           resdata{i}.(fn{1}) = resdata{i}.ModelStruct.Act.(fn{1});
        end
    end
       
    ClrMap     = [[31,119,180]./255;
                    [255,127,14]./255;
                    [31,119,180]./255;
                    [255,127,14]./255];
    % MarkerMap
    
    % LineStyleMap
    LSMap = ["-","-","--","--"];
    LW = [1.0,1.5];
    
    figID = 1;
    %% Control
    % Total torque commanded is similar
    % Total Control Action
    figID = figID+1;
    tVec = resdata{1}.tVec;
    NT = length(tVec);
    dt = tVec(2) - tVec(1);
    uNorm = zeros(N,NT);
    for i=1:N
        uNorm(i,:) = vecnorm(resdata{i}.UVec,1,2)';
    end
    
    tdiff = 150;
    IdConv = tVec > tdiff;
    unormnorm = sum(uNorm(:,IdConv).^2/length(IdConv),2).^0.5;
    
    printstr = 'Torque RMS: ';
    for i=1:N
        printstr = [printstr , ScenNames{i}, ': ',num2str(unormnorm(i)), ' Nm, '];
    end
    
    disp(printstr)
    
    unormmean = mean(uNorm(:,IdConv),2);
    printstr = 'Torque Mean l1: ';
    for i=1:N
        printstr = [printstr , ScenNames{i}, ': ',num2str(unormmean(i)), ' Nm, '];
    end
    disp(printstr)
    
    fig = figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    % Centralized
    hold on;
    title('Control Action','Interpreter','latex','fontsize',12)
    h = cell(1,N);
    for i=1:N
        h{i}=plot(tVec,vecnorm(resdata{i}.u_ref,1,1)','Color','k','LineWidth',LW(1));
        hold on;
    end
    h2 = cell(1,N);
    for i=1:N
        h2{i}=plot(tVec,vecnorm(resdata{i}.UVec,1,2)','Color',ClrMap(i,:),...
                                        'LineStyle', LSMap(i), 'LineWidth',LW(1));
        hold on;
    end
    hall = [h(end) h2(:)'];
    legend([hall{:}],['ref'; ScenNames(:)],'Location','northeast',...
                                            'Interpreter','latex','fontsize',12)
    ylabel('$\|\tau_k\|_1$ [Nm]','Interpreter','latex','fontsize',12)
    box on 
    xlabel('time [s]','Interpreter','latex','fontsize',12)
    savefig(fig,[OutputFile 'RWTNorm.fig']);
    saveas(fig,[OutputFile 'RWTNorm.png']);
    
    %% Control Integral
    figID = figID + 1;
    fig = figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    hold on;
    title('Control Action Integral','Interpreter','latex','fontsize',12)
    hold on;
    for i=1:N
        h2{i}=plot(tVec,cumsum(vecnorm(resdata{i}.UVec,1,2)'), 'Color',...
                           ClrMap(i,:),'LineStyle', LSMap(i), 'LineWidth',1.5);
        hold on;
    end
    hall = [h2(:)'];
    legend([hall{:}],[ScenNames(:)],'Location','northwest',...
                                            'Interpreter','latex','fontsize',12)
    ylabel('Total Torque [Nm]','Interpreter','latex','fontsize',12)
    box on 
    xlabel('time [s]','Interpreter','latex','fontsize',12)
    savefig(fig,[OutputFile 'RWTIntNorm.fig']);
    saveas(fig,[OutputFile 'RWTIntNorm.png']);
    
    
    %% LQR Cost Integral
    figID = figID + 1;
    fig = figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    % Centralized
    % Distributed
    h  = cell(1,N);
    h2 = cell(1,N);
    JInt   = zeros(N,NT);
    uRuInt = zeros(N,NT);
    for i=1:N
        QCtrl    = resdata{i}.SimParam.QCtrl;
        RCtrl    = resdata{i}.SimParam.RCtrl;
        T2Ctrl   = resdata{i}.OBSWStruct.Ctrl.T2Ctrl;
        TAct     = resdata{i}.OBSWStruct.Ctrl.Tact;
        Qsqt     = chol(QCtrl); 
        Rsqt     = chol(RCtrl);
        x_nom    = zeros(resdata{i}.IdStruct.RealMdl.NX,1);
        x_nom(resdata{i}.IdStruct.RealMdl.IdXRW) = resdata{i}.nomrwspeed;
        x_ref    = resdata{i}.x_ref;
        u_ref    = resdata{i}.u_ref;
        X        = resdata{i}.CLStates;
        UVec     = resdata{i}.UVec;
        xQx      = vecnorm(Qsqt * T2Ctrl * (X'-x_ref-x_nom),2,1).^2;
        uRu      = vecnorm(Rsqt * TAct * (UVec'-u_ref),2,1).^2;
        J        = 0.5 * (uRu + xQx);
        JInt(i,:)   = cumsum(J);
        uRuInt(i,:) = cumsum(0.5*uRu);
        h{i}   = semilogy(tVec,cumsum(J),'Color',ClrMap(i,:),...
                                       'LineStyle', '-','LineWidth',1.5);
        hold on;
        h2{i}  = semilogy(tVec,cumsum(0.5*uRu),'Color',[ClrMap(i,:), 0.5], ...
                                       'LineStyle', '--', 'LineWidth',1.5);
        hold on;
    end
    hall = [h(:)' h2(:)'];
    ScenNamesTC = strcat(ScenNames,' Total Cost');
    ScenNamesCC = strcat(ScenNames,' Control Cost');
    lgdall = [ScenNamesTC(:)' ScenNamesCC(:)'];
    
    legend([hall{:}],[lgdall(:)],'Location','northeast', 'Interpreter',...
                                                        'latex','fontsize',12)
    
    title('LQR Cost Function Integral','Interpreter','latex','fontsize',12)
    ylabel('LQR Cost Integral','Interpreter','latex','fontsize',12)
    box on 
    ylim([1,1e5])
    xlabel('time [s]','Interpreter','latex','fontsize',12)
    savefig(fig,[OutputFile 'LQRIntNorm.fig']);
    saveas(fig,[OutputFile 'LQRIntNorm.png']);
    
    % uRu J
    % Print final cost integral
    printstr  = 'LQR Total Cost Integral: ';
    printstr2 = 'LQR Control Cost Integral: ';
    for i=1:N
        printstr  = [printstr , ScenNames{i}, ': ',num2str(JInt(i,end)), ', '];
        printstr2 = [printstr2 , ScenNames{i}, ': ',num2str(uRuInt(i,end)), ', '];
    end
    
    disp(printstr)
    disp(printstr2)
    
    %% RW Speed
    figID = figID+1;
    % Centralized
    fig=figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    ylbls = {'$\delta \omega_{rw_x}$ [rpm]'};
    titles = ['Centralized X-axis $\delta \omega_{rw}$ wrt nominal';
              'Distributed X-axis $\delta \omega_{rw}$ wrt nominal'];
    for i=1:N
        subplot(1,N,i)
        ax = gca;
        ax.FontSize = 12;
        ax.TickLabelInterpreter='latex';
        minrwentmin = inf;
        maxrwentmax = -inf;
        IdS    = resdata{i}.IdStruct;
        Naxisc = min(6,max(IdS.RWperpoint));
    
        nomwrw = resdata{i}.nomrwspeed;
        maxwrw = resdata{i}.maxrwspeed;
        minwrw = resdata{i}.minrwspeed;
        maxrwent = maxwrw(1) - nomwrw(1);
        minrwent = minwrw(1) - nomwrw(1);
        X     = resdata{i}.CLStates;
        x_ref = resdata{i}.x_ref;
        rectangle(Position=[tVec(1),maxrwent*60/(2*pi), tVec(end)-tVec(1),...
                                           2*maxrwent*60/(2*pi)], ...
                                         FaceColor=ClrMap(i,:), ...
                                         FaceAlpha=.2,...
                                         EdgeColor=[ClrMap(i,:),.2])
        rectangle(Position=[tVec(1),2*minrwent*60/(2*pi), tVec(end)-tVec(1),...
                                               -minrwent*60/(2*pi)], ...
                                             FaceColor=ClrMap(i,:), ...
                                             FaceAlpha=.2,...
                                             EdgeColor=[ClrMap(i,:),.2])
        
        for j=1:3:Naxisc 
            hold on;
            for k=1:length(IdS.RWperpoint)
                Nrwm1 = sum(IdS.RWperpoint(1:k-1));
                if IdS.RWperpoint(k) >= j
                    wrw = X(:,IdS.RealMdl.IdXRW(Nrwm1+j))*60/(2*pi);
                    dwrw = wrw - nomwrw(Nrwm1+j)*60/(2*pi);
                    h =plot(tVec,x_ref(IdS.RealMdl.IdXRW(Nrwm1+j),:)*60/(2*pi),...
                        'Color','k','LineStyle',LSMap(i),'LineWidth',1.5);
                end
            end
        end
        for j=1:3:Naxisc 
            hold on;
            for k=1:length(IdS.RWperpoint)
                Nrwm1 = sum(IdS.RWperpoint(1:k-1));
                if IdS.RWperpoint(k) >= j
                    wrw = X(:,IdS.RealMdl.IdXRW(Nrwm1+j))*60/(2*pi);
                    dwrw = wrw - nomwrw(Nrwm1+j)*60/(2*pi);
                    h2=plot(tVec,dwrw','Color',ClrMap(i,:),'LineWidth',1.5,...
                                                        'LineStyle',LSMap(i));
                end
            end
            grid on;
            ylabel(ylbls{1},'interpreter','latex','fontsize',12)
            if j ==1
                title(titles(i,:),'Interpreter','latex','fontsize',12)
            end
        end
        minrwentmin = min([minrwentmin, minrwent]);
        maxrwentmax = max([maxrwentmax, maxrwent]);
        legend([h h2],{'$\delta\omega_{ref}$','$\delta\omega$'},...
                                'Location','northeast','Interpreter','latex','fontsize',12)
        xlabel('time [s]','fontsize',12)
        ylim(1.1*[minrwentmin,maxrwentmax]*60/(2*pi))
    end
    savefig(fig,[OutputFile 'WRWSep.fig']);
    saveas(fig,[OutputFile 'WRWSep.png']);
    
    %% RW Total Angular Momentum Along X axis - will be zero. Maybe norm?
    figID = figID+1;
    fig=figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    ylbls = {'X [rpm]','Y [rpm]','Z [rpm]'};
    h  = cell(1,N);
    h2 = cell(1,N);
    hold on;
    for i=1:N
        IdS   = resdata{i}.IdStruct;
        Naxisc   = min(6,max(IdS.RWperpoint));
        maxwrw   = resdata{i}.maxrwspeed;
        minwrw   = resdata{i}.minrwspeed;
        nomwrw   = resdata{i}.nomrwspeed;
        maxrwent = maxwrw(1) - nomwrw(1);
        minrwent = minwrw(1) - nomwrw(1);
        if i == 1
            rectangle(Position=[tVec(1),maxrwent*60/(2*pi), tVec(end)-tVec(1),...
                                                   2*maxrwent*60/(2*pi)], ...
                                                 FaceColor=[0,0,0], ...
                                                 FaceAlpha=.2,...
                                                 EdgeColor=[0,0,0,.2])
            rectangle(Position=[tVec(1),2*minrwent*60/(2*pi), tVec(end)-tVec(1),...
                                                   -minrwent*60/(2*pi)], ...
                                                 FaceColor=[0,0,0], ...
                                                 FaceAlpha=.2,...
                                                 EdgeColor=[0,0,0,.2])
        end
        ylbls = {'$\delta \omega_{rw}$ [rpm]'};
        X     = resdata{i}.CLStates; 
        x_ref = resdata{i}.x_ref;
        for j=1:3:Naxisc
            hold on;
            for k=1:length(IdS.RWperpoint)
                Nrwm1 = sum(IdS.RWperpoint(1:k-1));
                if IdS.RWperpoint(k) >= j
                    wrw = X(:,IdS.RealMdl.IdXRW(Nrwm1+j))*60/(2*pi);
                    dwrw = wrw - nomwrw(Nrwm1+j)*60/(2*pi);
                    h{i}  = plot(tVec,x_ref(IdS.RealMdl.IdXRW(Nrwm1+j),:)*60/(2*pi),'k','LineWidth',1.5);
                    h2{i} = plot(tVec,dwrw,'Color',ClrMap(i,:),'LineWidth',1.5);
                end
            end
        end
    end
    
    % grid on;
    box on
    ylabel(ylbls{1},'interpreter','latex','fontsize',12)
    title('X-axis $\delta \omega_{rw}$ wrt nominal','Interpreter','latex','fontsize',12)
    
    hall = [h(1)' h2(:)'];
    ScenNamesDW = strcat(ScenNames,' $\delta\omega_{rw}$');
    lgdall = ['$\delta\omega_{ref}$', ScenNamesDW(:)'];
    
    legend([hall{:}],[lgdall(:)],'Location','northeast', 'Interpreter',...
                                                        'latex','fontsize',12)
    xlabel('time [s]','interpreter','latex','fontsize',12)
    ylim(1.1*[minrwent,maxrwent]*60/(2*pi))
    savefig(fig,[OutputFile 'WRWAll.fig']);
    saveas(fig,[OutputFile 'WRWAll.png']);
    
    %% Attitude
    figID = figID+1;
    fig=figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    h = cell(1,N);
    h2 = cell(1,N);
    SetTime = zeros(1,N);
    printstr = 'Settling Time: ';
    for i=1:N
        IdS = resdata{i}.IdStruct;
        IdXrigrq = IdS.RealMdl.IdXrigrq;
        Nrigq    = IdS.RealMdl.Nrigrq;
        x_ref    = resdata{i}.x_ref;
        X        = resdata{i}.CLStates;
        AttX     = X(:,IdXrigrq)';
        AttX_ref = x_ref(IdXrigrq,:);
        h2{i} = plot(tVec, AttX(1,:)'*4/pi*180,'Color',ClrMap(i,:),...
                                          'LineStyle',LSMap(i),'LineWidth',1.5);
        hold on;
        % 
        Att_Tgt     = zeros(3,1);
        SlewRng     = resdata{i}.SimParam.TgtSlewAngleRng;
        Att_Tgt(1) =  SlewRng(1) + (SlewRng(2) - SlewRng(1));
        atterr = vecnorm(AttX*4/pi*180 - Att_Tgt,2,1);
        % compute energy along trajectory
        NX = resdata{i}.IdStruct.RealMdl.NX;
        IdXdrig  = IdS.RealMdl.IdXdrig;
        IdXd     = IdS.RealMdl.IdXd;
        Nfilt      = IdS.RealMdl.NFilt;
        IdXmodq  = IdS.RealMdl.IdXmodq;
        IdXRW    = IdS.RealMdl.IdXRW;
        IdXd = setdiff(IdXd,IdXRW);
        IdXJ  = [IdXmodq, IdXd];
        MassMat  = resdata{i}.ModelStruct.DescMat(IdXd-Nfilt,IdXd-Nfilt);  
        StiffMat     = resdata{i}.ModelStruct.StiffMatModes;
        H            = blkdiag(StiffMat,MassMat); % , diag(J_RW));
        nom_x = zeros(NX,1);
        nom_x(IdXRW) = resdata{i}.nomrwspeed;
        dX    = X(:,IdXJ)' - nom_x(IdXJ);
        JCL = zeros(1,NT);
        for j=1:NT
            JCL(j) = 0.5*dX(:,j)'*H*dX(:,j);
        end
        MRW = resdata{i}.RWDir.*resdata{i}.InertiaRW;
        MRW = MRW(:,1);
        nomrw = resdata{i}.nomrwspeed;
        maxrw = resdata{i}.maxrwspeed;
        dwrw = (maxrw- nomrw)*4/5;
        Jslew = 0.5 * sum(abs(MRW).*dwrw)^2 / MassMat(4,4);
        %Jslew = 0.5*sum(abs(MRW).*(dwrw.^2));
        Id    = find(and(JCL < Jslew * 0.02, atterr < Att_Tgt(1) * 0.02),1,'first');
        if Id
            SetTime(i) = tVec(Id)-resdata{i}.SimParam.TSlewStart;
            printstr = [printstr , ScenNames{i}, ': ',num2str(SetTime(i)), ' s, '];
        end
    end
    disp(printstr)
    
    box on
    hall = [h(end) h2(:)'];
    legend([hall{:}],['ref'; ScenNames(:)],'Location','southeast',...
                                            'Interpreter','latex','fontsize',12)
    ylabel('$$\theta_x$$ [deg]','Interpreter','latex','fontsize',12)
    title('Rotation Angle Around X','Interpreter','latex','fontsize',12)
    xlabel('time [s]','interpreter','latex','fontsize',12)
    savefig(fig,[OutputFile 'Attitude.fig']);
    saveas(fig,[OutputFile 'Attitude.png']);
    
    %% Attitude Error
    tdiff = 150; % 0.001; % 
    printstr = 'Pointing APE RMS: ';
    figID = figID+1;
    fig=figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    hold on;
    AttErr = zeros(N,NT);
    h = cell(1,N);
    for i=1:N
        IdXrigrq = resdata{i}.IdStruct.RealMdl.IdXrigrq;
        x_ref    = resdata{i}.x_ref;
        X        = resdata{i}.CLStates;
        quatref  = mrp2quat(x_ref(IdXrigrq,:)');
        quatref(~any(quatref~=0,2),1) = 1;
        z = quatmultiply(quatconj(quaternion(mrp2quat(X(:,IdXrigrq)))),...
                       quaternion(quatref)); 
        z = compact(z);
        AttErr(i,:) = abs(2*atan2d(vecnorm(z(:,2:4),2,2),z(:,1))*3600);
        IdConv = ceil(tdiff/(dt)):NT;
        aperms = sum(AttErr(i,IdConv).^2/length(IdConv)).^0.5;
        printstr = [printstr , ScenNames{i}, ': ',num2str(aperms), ' ", '];
        h{i} = plot(tVec, AttErr(i,:),'Color',ClrMap(i,:),...
                                     'LineStyle',LSMap(i),'LineWidth',1.5);
        hold on;
    end
    disp(printstr)
    
    box on
    title('Pointing Error wrt reference','Interpreter','latex','fontsize',12)
    ylabel('Point. Error [arcsec]','Interpreter','latex','fontsize',12)
    xlabel('time [s]','interpreter','latex','fontsize',12)
    
    legend([h{:}],[ScenNames(:)],'Location','northeast', 'Interpreter',...
                                                        'latex','fontsize',12)
    
    ylim([0, 1.1*max(AttErr(:),[],'all')]) 
    xlim([tVec(1),tVec(end)])
    % title('$$Attitude Quaternion$$','Interpreter','latex')
    savefig(fig,[OutputFile 'AttitudeErr.fig']);
    saveas(fig,[OutputFile 'AttitudeErr.png']);
    
    %% Angular Velocity Error?
    figID = figID+1;
    FirstMode = zeros(N,NT);
    
    tdiff = 150; % 0.01; % 150;
    IdConv = ceil(tdiff/(dt)):NT;
    errcl = zeros(1,N);
    printstr = 'RMS Ang Vel: ';
    fig=figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    h = cell(1,N);
    for i=1:N
        IdXrigrd  = resdata{i}.IdStruct.RealMdl.IdXrigrd;
        X = resdata{i}.CLStates;
        xref = resdata{i}.x_ref;
        FirstMode(i,:) = vecnorm(X(:,IdXrigrd)-xref(IdXrigrd,:)',2,2);
        errcl(i) = sum(FirstMode(IdConv).^2/length(IdConv)).^0.5;
        printstr = [printstr , ScenNames{i}, ': ',num2str(rad2deg(errcl(i))), ' deg/s, '];
        h{i}=plot(tVec, rad2deg(FirstMode(i,:)),'Color',ClrMap(i,:),'LineStyle',...
                                                    LSMap(i),'LineWidth',1.5);
        hold on;
    end
    disp(printstr)
    
    box on
    title('Angular Velocity Error wrt reference','Interpreter','latex','fontsize',12)
    ylabel('$$\|\omega\|$$ [deg/s]','Interpreter','latex','fontsize',12)
    xlabel('time [s]','Interpreter','latex','fontsize',12)
    legend([h{:}],[ScenNames(:)],'Interpreter','latex','fontsize',12)
    ylim([0, 1.1*max(rad2deg(FirstMode(:)),[],'all')])
    xlim([tVec(1),tVec(end)])
    savefig(fig,[OutputFile 'AngVelErr.fig']);
    saveas(fig,[OutputFile 'AngVelErr.png']);
    
    %% RMS Deformation
    figID = figID+1;
    
    fig=figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    
    tdiff = 150; % 0.01; % 
    errRMS = zeros(1,N);
    errPeak= zeros(1,N);
    printstr = 'RMS Deformation: ';
    printstr2 = 'Peak Deformation: ';
    
    J = zeros(N,NT);
    h = cell(1,N);
    for i=1:N
        IdXmodq    = resdata{i}.IdStruct.RealMdl.IdXmodq;
        X          = resdata{i}.CLStates;
        scmass     = resdata{i}.SCMass;
        J(i,:)     = vecnorm(X(:,IdXmodq),2,2)/(scmass^0.5);
        IdConv     = ceil(tdiff/dt):NT;
        errRMS(i)  = sum(J(i,IdConv).^2/length(IdConv)).^0.5;
        errPeak(i) = norm(J(i,:),'inf');
        printstr = [printstr , ScenNames{i}, ': ',num2str(1e6*errRMS(i)), ' microm, '];
        printstr2 = [printstr2 , ScenNames{i}, ': ',num2str(1e6*errPeak(i)), ' microm, '];
        h{i} = plot(tVec,1e6*J(i,:),'Color',ClrMap(i,:),'LineStyle',...
                                                    LSMap(i),'LineWidth',1.5);
        hold on
    end
    disp(printstr)
    disp(printstr2)
    box on
    ylim([0, 1.1e6*max(J(:))])
    ylabel('RMS Def. [$\mu m$]','Interpreter','latex','fontsize',12)
    legend([h{:}],[ScenNames(:)],'Interpreter','latex','fontsize',12)
    title('RMS Deformation [$\mu m$]','Interpreter','latex','fontsize',12)
    xlabel('time [s]','Interpreter','latex','fontsize',12)
    savefig(fig,[OutputFile 'RMSDefDistvsCent.fig']);
    saveas(fig,[OutputFile 'RMSDefDistvsCent.png']);
    
    
    %% Second Bending mode
    figID = figID+1;
    fig=figure(figID);
    ax = gca;
    ax.FontSize = 12;
    ax.TickLabelInterpreter='latex';
    hold on;
    FirstMode = zeros(N,NT);
    h = cell(1,N);
    for i=1:N
        IdXmodq   = resdata{i}.IdStruct.RealMdl.IdXmodq(2);
        FirstMode(i,:) = resdata{i}.CLStates(:,IdXmodq);
        h{i} = plot(tVec, FirstMode(i,:),'Color',ClrMap(i,:),'LineStyle',LSMap(i),...
                                                        'LineWidth',1.5);
        hold on;
    end
    % grid on;
    box on
    ylabel('$$\eta_2$$','Interpreter','latex','fontsize',12)
    xlabel('time [s]','Interpreter','latex','fontsize',12)
    legend([h{:}],[ScenNames(:)],...
                        'Interpreter','latex','fontsize',12)
    title('Second Bending Mode','Interpreter','latex','fontsize',12)
    ylim([1.1*min(FirstMode(:),[],'all'), 1.1*max(FirstMode(:),[],'all')])
    xlim([tVec(1),tVec(end)])
    savefig(fig,[OutputFile 'scndbendmode.fig']);
    saveas(fig,[OutputFile 'scndbendmode.png']);

end

