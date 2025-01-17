function PlotCampbellDiagram(figID, sysMeas, ModelStruct, SimParam, IdStruct, ...
                                                    wrwnomVec, SVRW, file)
    % compute singular values for plot (log scale to look nice, sparser to
    % not take too long
    NrwPcd       = 100;
    Nfreq        = 2e3; % 1e4; % 
    freqmin      = SimParam.TgtMinFreq/(2*pi); 
    % (SimParam.TgtMinBdwt/2)/(2*pi); % 1e-3; % 0.6; % 0.1;
    freqmax      = SimParam.SamplingFreq/2;
    plotfreqcdmin  = (SimParam.TgtMinBdwt/2)/(2*pi); 
    plotfreqcdmax  = SimParam.SamplingFreq/2;
    wrwnomcdVec    = linspace(freqmin,freqmax,NrwPcd);
    wsigmacd       = logspace(log10(plotfreqcdmin), log10(plotfreqcdmax), Nfreq)*(2*pi);
    NatFreqs     = zeros(IdStruct.RealMdl.NX,NrwPcd);
    inst         = zeros(NrwPcd,     1);
    SVs          = zeros(Nfreq, NrwPcd);
    IdUWRWW      = IdStruct.IdUWRWW;
    Cnew         = eye(IdStruct.RealMdl.NX);
    Cnew         = Cnew(IdStruct.RealMdl.IdXmodd,:);
    Dnew         = zeros(IdStruct.RealMdl.Nmod, IdStruct.NUW);
    RWDir        = ModelStruct.Act.RWDir;
    InertiaRW    = ModelStruct.Act.InertiaRW;
    Nrw          = ModelStruct.Act.Nrw;
    DampMatModes = ModelStruct.DampMatModes;
    WPert       = ModelStruct.Pert.WPert;
    % Campbell diagram plot 
    sysMeas0 = dss(sysMeas.A,sysMeas.B,Cnew,Dnew,sysMeas.E);
    for irw = 1:NrwPcd
        wrw = wrwnomcdVec(irw)*2*pi;
        % RW Gyroscopic torque
        Jrwwrw = zeros(3,IdStruct.RealMdl.Nrigrd + IdStruct.RealMdl.Nmod);
        Nrwrw  = zeros(IdStruct.RealMdl.Nmod, IdStruct.RealMdl.Nrigrd ...
                                                   + IdStruct.RealMdl.Nmod);
        for iu = 1:IdStruct.NU
            jrw0 = skew(RWDir(iu,:)' * (InertiaRW(iu) .* wrw));
            Jrwwrw = Jrwwrw + jrw0 * [eye(3), squeeze(Nrw(iu,:,:))];
            Nrwrw  = Nrwrw  + squeeze(Nrw(iu,:,:))' *  jrw0 * [eye(3), squeeze(Nrw(iu,:,:))];
        end
        sysMeas0.A(IdStruct.RealMdl.IdXrigrd, IdStruct.RealMdl.IdXrigrd) = ...
                                                           Jrwwrw(:,1:3);
        % adding gyricity
        sysMeas0.A(IdStruct.RealMdl.IdXrigrd, IdStruct.RealMdl.IdXmodd) = ...
                                                           Jrwwrw(:,4:end);
        sysMeas0.A(IdStruct.RealMdl.IdXmodd, IdStruct.RealMdl.IdXrigrd) = ...
                                          Nrwrw(:,1:3);
        sysMeas0.A(IdStruct.RealMdl.IdXmodd, IdStruct.RealMdl.IdXmodd) = ...
                                          -DampMatModes + Nrwrw(:,4:end);

        NPertperNU  = 2*SimParam.RadFPertRW + 2*SimParam.RadMPertRW;
        idW = NPertperNU*(iu-1)+(1:NPertperNU); 
        SVs(:,irw) = max(sigma(sysMeas0(:,IdUWRWW(idW))*WPert(idW,idW),wsigmacd),[],1)';
        [Wn,~] = damp(sysMeas0.E\sysMeas0.A);
        inst(irw) = max(real(eig(sysMeas0)));
        NatFreqs(:,irw) = sort(Wn)/(2*pi);
    end
    
    x2 = [wrwnomVec-(SimParam.TgtMinBdwt/(4*pi)); ...
         wrwnomVec+(SimParam.TgtMinBdwt/(4*pi))]*60;
    curve1 = wrwnomVec; % -(SimParam.TgtMinBdwt/2)/(2*pi);
    curve2 = wrwnomVec; % +(SimParam.TgtMinBdwt/2)/(2*pi);
    inBetween = [curve1; curve2];
    cmap2 = flipud(hot(10000));
    % 
    % Plot the line with width 8 so we can see the colors well.
    
    % Campbell diagram
    fig=figure(figID);
    % Sigma values as a function of RW Speed
    ax1 = axes;
    ax1.TickLabelInterpreter='latex';
    ax1.FontSize = 12;
    
    ha=contourf(wrwnomcdVec*60,wsigmacd/(2*pi),mag2db(SVs),500,'EdgeAlpha',0);
    hold on;
    % colormap(flipud(gray));
    colormap(flipud(bone));
    cmin = min(mag2db(SVs),[],'all');
    cmax = max(mag2db(SVs),[],'all');
    clim([cmin, cmin + (cmax - cmin)*1.5])
    hc=loglog(wrwnomcdVec*60,NatFreqs,'k');
    hold on;
    % he=plot(ModelStruct.Act.nomrwspeed(iu)*60/(2*pi), ...
    %             ModelStruct.Act.nomrwspeed(iu)/(2*pi),'bx','LineWidth',1.5);


    % change color to match convolution of svd
    NrwSV = length(wrwnomVec);
    ax2 = axes;
    ax1.TickLabelInterpreter='latex';
    ax1.FontSize = 12;
    hb=surface(x2, inBetween, ...
        [zeros(2,NrwSV)], mag2db([SVRW;SVRW]),...
    'EdgeColor', 'interp','LineWidth', 1.25,'EdgeAlpha',0.15,'FaceAlpha',0.15);% 'FaceColor', 'no', , 'EdgeAlpha', 'interp'
    % linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    hold on;

    % plot bar and marker on selected nominal velocity range
    
    % hold on;
    % hd = fill(x2, inBetween, 'b','EdgeColor','none','FaceAlpha',0.3);
    set(ax1, 'XScale', 'log', 'YScale', 'log', 'XScale', 'log');
    xlabel(ax1,'Wheel Speed [rpm]','interpreter','latex','fontsize',12)
    ylabel(ax1,'Frequency [Hz]','interpreter','latex','fontsize',12)
    ylim(ax1,[(SimParam.TgtMinBdwt/2)/(2*pi), ...
        wsigmacd(end)./(2*pi)])
    xlim(ax1,[wrwnomcdVec(1), wrwnomcdVec(end)]*60)
    cb = colorbar(ax1);
    cb.Label.String = 'System Response to RW Perturbations [dB]';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter='latex';
    cb.Label.FontSize = 12;
    cb.FontSize    = 12;
    axis(ax2,[[wrwnomcdVec(1), wrwnomcdVec(end)]*60,...
                                [wsigmacd(1), wsigmacd(end)]./(2*pi)])
    set(ax2, 'XScale', 'log', 'YScale', 'log', 'XScale', 'log');
    colormap(ax2, cmap2)
    clim(ax2,mag2db([min(SVRW), max(SVRW)]))
    c2 = colorbar(ax2,'northoutside', 'AxisLocation','out');
    c2.Label.String = 'Main harmonic system response magnitude [dB]';
    c2.Label.Interpreter = 'latex';
    c2.TickLabelInterpreter='latex';
    c2.Label.FontSize    = 12;
    c2.FontSize    = 12;
    linkprop([ax1, ax2], {'View', 'PlotBoxAspectRatio', 'Position', ...
        'XLim', 'YLim','YScale'});
    xticks(ax1,[200,300,500,1000,2000,4000,6000])
    legend([hc(1)],{'Natural Frequencies'}, ...
                        'location','southeast','Interpreter','latex','fontsize',12)
    % legend([hc(1), he],{'Natural Frequencies','Optimal RW Nom. Speed'}, ...
    %                     'location','southeast','Interpreter','latex','fontsize',12)
    savefig(fig,[file 'CampDiag.fig']);
    saveas(fig,[file 'CampDiag.png']);


end
