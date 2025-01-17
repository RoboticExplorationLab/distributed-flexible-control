function [TKD, IdKD]    = AdHocKD(SysIddModel, ModelStruct, IdStruct, SimParam)
    %ADHOCKD Summary of this function goes here
    %   Detailed explanation goes here
    Mrw = (ModelStruct.Act.InertiaRW .* ModelStruct.Act.RWDir);
    
    if IdStruct.NU > 3
        tctrl    = Mrw./repmat(vecnorm(Mrw),IdStruct.NU,1);
        RWCtrMat = [tctrl';null(tctrl')'];
        T2ctrl   = eye(IdStruct.RealMdl.NX);
        T2ctrl(IdStruct.RealMdl.IdXRW, IdStruct.RealMdl.IdXRW) = RWCtrMat;
    else
        RWCtrMat = eye(IdStruct.RealMdl.NXRW);
        T2ctrl   = eye(IdStruct.RealMdl.NX);
    end
    IdXrigt  = [IdStruct.RealMdl.IdXrigtq, IdStruct.RealMdl.IdXrigtd];
    % T2ctrl(IdXrigt,:) = [];
    SysIddModelCtrOnly = ss2ss(SysIddModel,T2ctrl);
    SysIddModelCtrOnly = SysIddModelCtrOnly * RWCtrMat;
    if IdStruct.NU > 3
        ObsUnctrlX         = [IdStruct.RealMdl.IdXRW(1:3)];
    else
        ObsUnctrlX = [];
    end
    UnobUnctUX         = [IdStruct.RealMdl.IdXrigtq, ...
                          IdStruct.RealMdl.IdXrigtd];
    UnctrlX            = [UnobUnctUX, ...
                          ObsUnctrlX];
    
    CaO    = T2ctrl(setdiff(1:end,UnctrlX),:)';
    UbO    = T2ctrl(ObsUnctrlX,:)';
    CbU    = [];
    UaU    = T2ctrl(UnobUnctUX,:)';
    TKD    = [CaO UbO CbU UaU];

    % Indexing of Kalman Decomposition
    % obs/unobs, ctrl/unctrl
    IdKD         = IdStruct.RealMdl;
    IdKD.NCOaOB  = size(CaO, 2);
    IdKD.NUCbOB  = size(UbO, 2);
    IdKD.NCObUO  = size(CbU, 2);
    IdKD.NUCaUO  = size(UaU, 2);
    IdKD.NOB     = IdKD.NCOaOB + IdKD.NUCbOB;
    IdKD.NUO     = IdKD.NCObUO + IdKD.NUCaUO;
    IdKD.NCO     = IdKD.NCOaOB + IdKD.NCObUO;
    IdKD.NUC     = IdKD.NUCbOB + IdKD.NUCaUO;
    IdKD.IdCOaOB = 1:IdKD.NCOaOB;
    IdKD.IdUCbOB = IdKD.NCOaOB + (1:IdKD.NUCbOB);
    IdKD.IdCObUO = IdKD.NCOaOB + IdKD.NUCbOB + (1:IdKD.NCObUO);
    IdKD.IdUCaUO = IdKD.NCOaOB + IdKD.NUCbOB + IdKD.NCObUO + (1:IdKD.NUCaUO);
    IdKD.IdOB    = [IdKD.IdCOaOB, IdKD.IdUCbOB];
    IdKD.IdUO    = [IdKD.IdCObUO, IdKD.IdUCaUO];
    IdKD.IdCO    = [IdKD.IdCOaOB, IdKD.IdCObUO];
    IdKD.IdUC    = [IdKD.IdUCbOB, IdKD.IdUCaUO];
    % actual state allocation
    % total number of position states and indices within them
    IdKD.IdXrigtq = IdKD.IdUCaUO(1:IdKD.Nrigtq);
    IdKD.IdXrigrq = IdKD.NFilt +  (1:IdKD.Nrigrq);
    IdKD.IdXrigq  = [IdKD.IdXrigtq, IdKD.IdXrigrq];
    IdKD.IdXmodq  = IdKD.NFilt + IdKD.Nrigrq + (1:IdKD.Nmod);
    IdKD.IdXq     = [IdKD.IdXrigq, IdKD.IdXmodq];
   
    % total number of velocity states and indices within them
    IdKD.NMomRW   = 3;
    IdKD.NNullRWd = IdKD.NXRW - IdKD.NMomRW;

    IdKD.IdXrigtd = IdKD.IdUCaUO(IdKD.Nrigtq+1:end);
    IdKD.IdXrigrd = IdKD.NFilt + IdKD.Nrigrq + IdKD.Nmod + (1:IdKD.Nrigrd);
    IdKD.IdXmodd  = IdKD.NFilt + IdKD.Nrigrq + IdKD.Nmod + IdKD.Nrigrd + (1:IdKD.Nmod);
    IdKD.IdXdrigr = [IdKD.IdXrigrq, IdKD.IdXrigrd];
    IdKD.IdXNullRWd = IdKD.NFilt + IdKD.Nrigrq + 2*IdKD.Nmod + IdKD.Nrigrd + (1:IdKD.NNullRWd);
    IdKD.IdXrigd  = [IdKD.IdXrigtd, IdKD.IdXrigrd];
    IdKD.IdXd     = [IdKD.IdXrigd, IdKD.IdXmodd, IdKD.IdXNullRWd];
    
    % unctrl
    IdKD.NMomRW   = 3;
    IdKD.IdXMomRW = IdKD.IdUCbOB;
    IdKD.IdXRW    = [IdKD.IdXNullRWd, IdKD.IdXMomRW];

    % unobs and unctrl
    IdKD.IdXdrigt = [IdKD.IdXrigtq, IdKD.IdXrigtd];
    % total number and Id of rigid body coordinates
    IdKD.IdXmod    = [IdKD.IdXmodq, IdKD.IdXmodd];
    IdKD.NX        = IdKD.NAAFilt + IdKD.NXq + IdKD.NXd;
end

