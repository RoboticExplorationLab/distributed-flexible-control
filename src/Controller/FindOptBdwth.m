function [ModelStruct, sysMeas] = FindOptBdwth(figID, sysMeas0, ModelStruct, ...
                                                  SimParam, IdStruct, file)

    RWDir        = ModelStruct.Act.RWDir;
    InertiaRW    = ModelStruct.Act.InertiaRW;
    Nrw          = ModelStruct.Act.Nrw;
    WPert        = ModelStruct.Pert.WPert;
    StateMat     = ModelStruct.StateMat;
    DescMat      = ModelStruct.DescMat;
    DampMatModes = ModelStruct.DampMatModes;
    IdUWRWW      = IdStruct.IdUWRWW;
    IdXmodd      = IdStruct.RealMdl.IdXmodd;
    NU           = IdStruct.NU;
    RealId       = IdStruct.RealMdl;
    NPertperNU  = 2*SimParam.RadFPertRW + 2*SimParam.RadMPertRW;

    % compute the singular values along the RW in a refined manner
    % frequencies need to match
    NrwP        = 1e3; % 1e4; % 
    freqmin     = SimParam.TgtMinFreq/(2*pi);
    freqmax     = SimParam.SamplingFreq/2;
    wrwnomVec   = linspace(freqmin,freqmax,NrwP);
    SVRW        = zeros(IdStruct.NU,NrwP);
    B           = sysMeas0.B;

    Cnew        = eye(IdStruct.RealMdl.NX);
    Cnew        = Cnew(IdStruct.RealMdl.IdXmodd,:);
    Dnew        = zeros(IdStruct.RealMdl.Nmod, IdStruct.NUW);
    sysMeas     = dss(sysMeas0.A,sysMeas0.B,Cnew,Dnew,sysMeas0.E);
    dwrw        = wrwnomVec(2)-wrwnomVec(1);
    NS          = round(SimParam.TgtMinBdwt/(2*pi*dwrw));
    NM          = ceil(NS/2);
    nomrwhz     = zeros(IdStruct.NU,1);
    nomrwrd     = zeros(IdStruct.NU,1);
    WPert       = ModelStruct.Pert.WPert;
    for irw = [1:IdStruct.NU,1]
        idW = NPertperNU*(irw-1)+(1:NPertperNU); 
        % needs to multiply with wrw
        SVRW(irw,:)  = max(sigma(sysMeas(:,IdUWRWW(idW))*WPert(idW,idW),wrwnomVec*2*pi),[],1);
        SVRW(irw,:)  = SVRW(irw,:) .* (wrwnomVec.^2);
        % MeanSVRW    = movmean(SVRW,NS);
        MeanSVRW     = movmax(SVRW(irw,:),NS);
        MeanSVRW2    = MeanSVRW(NM:end-NM);
        wrwnomVec2   = wrwnomVec(NM:end-NM);
        [~,Idxmean]  = min(MeanSVRW2);%  + NM - 1;
        nomrwhz(irw) = wrwnomVec2(Idxmean);
        nomrwrd(irw) = wrwnomVec2(Idxmean)*(2*pi);
    end
    maxrw       = nomrwhz+(SimParam.TgtMinBdwt/2)/(2*pi);
    minrw       = nomrwhz-(SimParam.TgtMinBdwt/2)/(2*pi);
    
    % Update the values
    % RW Gyroscopic torque
    ModelStruct.Act.nomrwspeed = nomrwrd;
    ModelStruct.Act.maxrwspeed = maxrw*2*pi;
    ModelStruct.Act.minrwspeed = minrw*2*pi;

    Jrwwrw = zeros(3,RealId.Nrigrd + RealId.Nmod);
    Nrwrw  = zeros(RealId.Nmod, RealId.Nrigrd + RealId.Nmod);
    for iu = 1:IdStruct.NU
        jrw0 = skew(RWDir(iu,:)' * (InertiaRW(iu) .* nomrwrd(iu)));
        omegarw = squeeze(Nrw(iu,:,:));
        omegarw = [eye(3), omegarw];
        Jrwwrw = Jrwwrw + jrw0 * omegarw;
        Nrwrw  = Nrwrw  + squeeze(Nrw(iu,:,:))' *  jrw0 * omegarw;
    end

    % adding gyricity
    sysMeas0.A(RealId.IdXrigrd, RealId.IdXrigrd) = Jrwwrw(:,1:3);
    sysMeas0.A(RealId.IdXrigrd, RealId.IdXmodd) = Jrwwrw(:,4:end);
    sysMeas0.A(RealId.IdXmodd, RealId.IdXrigrd) = Nrwrw(:,1:3);

    sysMeas0.A(RealId.IdXmodd, RealId.IdXmodd) = ...
                                      -DampMatModes + Nrwrw(:,4:end);
    sysMeas = sysMeas0;
    figID = figID + 1;
    PlotCampbellDiagram(figID, sysMeas, ModelStruct, SimParam,...
                                IdStruct, wrwnomVec, MeanSVRW,file)
    
    % first RW bandwidth optimization
    % turn into separate function
    % make a plot for individual RW speed optimization 
    figID                     = figID + 1;
    PlotConfig.plotfreqsvdmin = 1e-1; % 1e-3; 
    PlotConfig.plotfreqsvdmax = 1e3;
    PlotConfig.NSig           = 1e3; % 1e4;
    PlotConfig.freqtypacs     = 5;
    PlotConfig.ftsize         = 12;
    PlotConfig.freqacs        = SimParam.SamplingFreq/2;
    PlotConfig.nomrw          = ModelStruct.Act.nomrwspeed/(2*pi);
    PlotConfig.minrw          = ModelStruct.Act.minrwspeed/(2*pi);
    PlotConfig.maxrw          = ModelStruct.Act.maxrwspeed/(2*pi);

    Cnew        = eye(IdStruct.RealMdl.NX);
    Cnew        = Cnew(IdStruct.RealMdl.IdXmodd,:);
    Dnew        = zeros(IdStruct.RealMdl.Nmod, IdStruct.NUW);
    sysMeasPlot = dss(sysMeas.A, sysMeas.B(:,IdUWRWW), Cnew, ...
                                            Dnew(:,IdUWRWW), sysMeas.E);
    
    PlotSingleOptimRWBdwth(sysMeasPlot(:,1:NPertperNU), PlotConfig, file, figID);
    figID = figID + 1;
    PlotAllOptimRWBdwth(sysMeasPlot,PlotConfig, IdStruct, file, figID);
end
