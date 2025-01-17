function PlotAllOptimRWBdwth(sysMeas,PlotConfig, IdStruct, file, figID)
    plotfreqsvdmin = PlotConfig.plotfreqsvdmin;
    plotfreqsvdmax = PlotConfig.plotfreqsvdmax;
    NSig           = PlotConfig.NSig;
    freqtypacs     = PlotConfig.freqtypacs;
    ftsize         = PlotConfig.ftsize;
    freqacs        = PlotConfig.freqacs;
    nomrw          = PlotConfig.nomrw;
    minrw          = PlotConfig.minrw;
    maxrw          = PlotConfig.maxrw;

    wsigmaplot = logspace(log10(plotfreqsvdmin),log10(plotfreqsvdmax),NSig)*(2*pi);
    SVplot = sigma(sysMeas,wsigmaplot);
    maxSVplot = max(SVplot,[],1);
    SVs = zeros(IdStruct.NU, length(wsigmaplot));
    for i=1:IdStruct.NU
        IdUWRWW = (1+(i-1)*4):(4*i);
        SVs(i,:) = max(sigma(sysMeas(:,IdUWRWW),wsigmaplot),[],1);
    end
    % max
    fig=figure(figID);
    semilogx(wsigmaplot/(2*pi), NaN(size(wsigmaplot)));
    hold on;
    ax=gca;
    ax.FontSize=ftsize;
    ax.TickLabelInterpreter='latex';
    cmin = min(mag2db(SVs),[],'all');
    cmax = max(mag2db(SVs),[],'all');
    hold on
    grayscale = [linspace(1,0,256)',linspace(1,0,256)',linspace(1,0,256)'];
    for i=1:IdStruct.NU
        % Line plots
        y = repmat((i-1:0.1:i)',1,length(wsigmaplot));
        x = repmat(wsigmaplot/(2*pi),size(y,1),1);
        SVu = repmat(SVs(i,:),size(y,1),1);
        ha=contourf(gca,x,y,mag2db(SVu),500, 'EdgeAlpha',0);
        hold(gca,'on');
        
        % RW bandwidth
        hb=patch(gca,'XData',[minrw(i);minrw(i);maxrw(i);maxrw(i)],'YData',[i-1;i;i;i-1], ...
            FaceColor=[0 0 0.9], FaceAlpha=.1, EdgeColor=[0 0 0.9]);

        % RW nominal speed
        hf = plot([nomrw(i), nomrw(i)],[i-1,i],'--','Color',[0,0,0.9]);
        % h3=xline(nomrw,'b--');
        % xline(maxrw(1),'b')
        % h5=xline(minrw,'b',{'RW Bandwidth'},'LabelHorizontalAlignment','left');
        % h5.FontSize = ftsize;
        % horizontal line below
        yline(i,'k')
    end
    colormap(gca,grayscale);
    clim(gca,[cmin, cmin + (cmax - cmin)*1.3])
    ylim(gca,[0, IdStruct.NU])
    % xlim(gca,[1e-1/(2*pi), freqacs])
    xlim(gca,[wsigmaplot(1)/(2*pi), wsigmaplot(end)/(2*pi)])
    ylbls = cell(1,IdStruct.NU);
    for i=1:IdStruct.NU
        if floor((i-1)/3) == 0
            dir = '+';
        else
            dir = '-';
        end
        if mod(i,3) == 1
            nb = strlength([num2str(i+1) ' Y']);
            ylbls{i} = ['RW ' dir 'X'];
        elseif mod(i,3) == 2
            ylbls{i} = ['RW ' dir 'Y'];
        else
            nb = strlength([num2str(i-1) dir 'Y']);
            ylbls{i} = ['RW ' dir 'Z'];
        end
    end
    set(gca,'YTick',(1:IdStruct.NU)-0.5,'YTickLabel',ylbls, 'XScale', 'log','YTickLabelRotation',45);
    % ADCS controlled 
    hold on;
    hc=patch(gca,'XData',[plotfreqsvdmin;plotfreqsvdmin;freqtypacs;...
        freqtypacs], 'YData',[0;IdStruct.NU;IdStruct.NU;0], ...
            FaceColor=[0.9 0.9 0.9], FaceAlpha=.1, EdgeColor=[0 0 0]);
    % xregion(plotfreqsvdmin,freqtypacs,FaceColor=[0.1 0.1 0.1 0.3],FaceAlpha=.3,...
    %                                                 EdgeColor=[0 0 0 1]);
    % extended ADCS controlled
    % hold on;
    hd=patch(gca,'XData',[freqtypacs;freqtypacs;freqacs;...
        freqacs], 'YData',[0;IdStruct.NU;IdStruct.NU;0], ...
            FaceColor=[0.9 0.9 0], FaceAlpha=.1, EdgeColor=[0.9 0.9 0]);
    % xregion(freqtypacs,freqacs,FaceColor=[0.1 0.1 0.1 0.1],FaceAlpha=.1,...
    %                                                 EdgeColor=[0 0 0 1]);
    % uncontrolled
    hold on;
    he=patch('XData',[freqacs;freqacs;plotfreqsvdmax;...
        plotfreqsvdmax], 'YData',[0;IdStruct.NU;IdStruct.NU;0], ...
            FaceColor=[0.9 0 0], FaceAlpha=.1, EdgeColor=[0.9 0 0]);
    % xregion(freqacs,plotfreqsvdmax,FaceColor=[0.9 0 0, 0.1],FaceAlpha=.1,...
    %                                                 EdgeColor=[0.4 0 0.7]);
    % ACS Controller Bandwidth
    h2=xline(freqacs,'k','LineWidth',1.5);%,{'Extended ACS Bandwidth/','Active Vibration Suppression'});
    h2.FontSize = ftsize;
    % ACS Controller Bandwidth
    h4=xline(freqtypacs,'k','LineWidth',1.5);%,{'ACS Bandwidth'});
    h4.FontSize = ftsize;
    h6=xline(plotfreqsvdmax-1e-6,'r','LineWidth',1.5);%,{'Uncontrolled'},'LabelHorizontalAlignment','left');
    h6.FontSize = ftsize;
    xlabel('Frequency [Hz]','interpreter','latex','FontSize',ftsize);
    % ylabel('RW','interpreter','latex','FontSize',ftsize);
    legend([hc hd he hb hf],{'ACS Bandwidth', sprintf('Extended ACS/\nVS Bandwidth'), ...
            'Uncontrolled Bdwth', 'RW Speed Constraint', 'Nominal RW Speed'}, ...
            'Location','southeast','interpreter','latex','FontSize',ftsize);
    % legend([h1 h3],{'Max. Singular Value', 'Nominal RW Speed'},...
    %                                              'Location','southeast')
    %  h2'RW Int. MA',
    savefig(fig,[file 'PSDSYSRWAll.fig']);
    saveas(fig,[file 'PSDSYSRWAll.png']);
end

