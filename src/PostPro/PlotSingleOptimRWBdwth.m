function PlotSingleOptimRWBdwth(sysMeas,PlotConfig, file, figID)
    % Currently picking the same bandwidth for all
    % max SV or max out of Bode plots for single input?
    % wsigma = linspace(10^(-3),100,1e5)*(2*pi);
    plotfreqsvdmin = PlotConfig.plotfreqsvdmin;
    plotfreqsvdmax = PlotConfig.plotfreqsvdmax;
    NSig           = PlotConfig.NSig;
    freqtypacs     = PlotConfig.freqtypacs;
    ftsize         = PlotConfig.ftsize;
    freqacs        = PlotConfig.freqacs;
    nomrw          = PlotConfig.nomrw(1);
    minrw          = PlotConfig.minrw(1);
    maxrw          = PlotConfig.maxrw(1);


    wsigmaplot = logspace(log10(plotfreqsvdmin),log10(plotfreqsvdmax),NSig)*(2*pi);
    SVplot = sigma(sysMeas,wsigmaplot);
    maxSVplot = max(SVplot,[],1);
    % max
    fig=figure(figID);
    semilogx(wsigmaplot/(2*pi), NaN(size(wsigmaplot)));
    hold on;
    ax=gca;
    ax.FontSize=ftsize;
    ax.TickLabelInterpreter='latex';
    % Line plots
    % y/u(s)=G(s)
    h1=semilogx(wsigmaplot/(2*pi), mag2db(maxSVplot),'k','LineWidth',0.75);
    % ADCS controlled
    xregion(plotfreqsvdmin,freqtypacs,FaceColor=[0.1 0.1 0.1 0.3],FaceAlpha=.3,...
                                                    EdgeColor=[0 0 0 1]);
    % extended ADCS controlled
    xregion(freqtypacs,freqacs,FaceColor=[0.9 0.9 0 0.1],FaceAlpha=.1,...
                                                    EdgeColor=[0 0 0 1]);
    % uncontrolled
    xregion(freqacs,plotfreqsvdmax,FaceColor=[0.9 0 0, 0.1],FaceAlpha=.1,...
                                                    EdgeColor=[0.4 0 0.7]);
    % RW bandwidth
    xregion(minrw,maxrw,FaceColor=[0 0 0.9 0.1],FaceAlpha=.1,...
                                                    EdgeColor=[0 0 0.9 1]);
    
    % h2=semilogx(wsigma(NM:end-NM)/(2*pi), mag2db(M),'b');
    % moving average of cost over single spectrum (ignoring gyricity)
    % semilogx(nomrw,mag2db(min(M)),'bo','LineWidth',1.5)
    % ACS Controller Bandwidth
    h2=xline(freqacs,'k',{'Extended ACS Bandwidth/','Active Vibration Suppression'});
    h2.FontSize = ftsize;
    % ACS Controller Bandwidth
    h4=xline(freqtypacs,'k',{'ACS Bandwidth'},'LabelVerticalAlignment','bottom');
    h4.FontSize = ftsize;
    % RW nominal speed
    h3=xline(nomrw,'b--');
    xline(maxrw(1),'b')
    h5=xline(minrw,'b',{'RW Speed Constraint'},'LabelHorizontalAlignment',...
                                'left','LabelVerticalAlignment','bottom');
    h5.FontSize = ftsize;
    h6=xline(plotfreqsvdmax,'r',{'Uncontrolled'},'LabelHorizontalAlignment','left');
    h6.FontSize = ftsize;
    xlabel('Frequency [Hz]','Interpreter','latex','FontSize',ftsize);
    ylabel('Magnitude [dB]','Interpreter','latex','FontSize',ftsize);
    legend([h1 h3],{'$\|G\|_{\infty}$', 'Nominal RW Speed'},...
                   'Location','southeast','Interpreter','latex','FontSize',ftsize);
    %  h2'RW Int. MA',
    savefig(fig,[file 'PSDSYSRW.fig']);
    saveas(fig,[file 'PSDSYSRW.png']);
end

