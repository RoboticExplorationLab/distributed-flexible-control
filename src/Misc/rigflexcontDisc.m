function xdot = rigflexcontDisc(xtrw,u, RealStruct, IdStruct, ModelStruct, SimParam)
    % load parameters
    m      = RealStruct.m;
    J      = RealStruct.J;
    Jrw    = RealStruct.Jrw;
    Leta   = RealStruct.Leta;
    Lr     = RealStruct.Lr;
    Minv   = RealStruct.Minv;
    MinvB  = RealStruct.MinvB;
    MinvK  = RealStruct.MinvK;
    MinvC  = RealStruct.MinvC;
    MinvG  = RealStruct.MinvG;
    EDir   = RealStruct.EDir;
    Nrw    = RealStruct.Nrw;
    FiltFlag = RealStruct.FiltFlag;
    % load indexes
    NXRW     = IdStruct.RealMdl.NXRW;
    Nrigrd   = IdStruct.RealMdl.Nrigrd;
    Nrigtd   = IdStruct.RealMdl.Nrigtd;
    Nmod     = IdStruct.RealMdl.Nmod;
    NX       = IdStruct.RealMdl.NX;
    NXd      = IdStruct.RealMdl.NXd;
    IdXrigrq = IdStruct.RealMdl.IdXrigrq;
    IdXrigtd = IdStruct.RealMdl.IdXrigtd;
    IdXrigrd = IdStruct.RealMdl.IdXrigrd;
    IdXdrigr = IdStruct.RealMdl.IdXdrigr;
    IdXRW    = IdStruct.RealMdl.IdXRW;
    IdXmodq  = IdStruct.RealMdl.IdXmodq;
    IdXmodd  = IdStruct.RealMdl.IdXmodd;
    IdXdmod  = IdStruct.RealMdl.IdXdmod;
    IdXFilt  = IdStruct.RealMdl.IdXFilt;
    % load states
    x       = xtrw(1:NX);
    thetarw = xtrw(NX+1:end);
    xdot    = zeros(NX+NXRW,1);
    p       = x(IdXrigrq);
    vsc     = x(IdXrigtd);
    wsc     = x(IdXrigrd);
    wrw     = x(IdXRW);
    eta     = x(IdXmodq);
    etadot  = x(IdXmodd);
    pp      = 1 + (p'*p);
    sp      = skew(p);
    % Linear momentum
    Mv        = m * vsc + Lr * etadot;
    % Angular momentum
    Mw        = J * wsc + Leta * etadot; %  + EDir' * (Jrw .* wrw);
    Mwrw      = (Jrw .* wrw) .* EDir;
    % Coriolis and Centrifugal Forces
    skewwsc = skew(wsc);
    Cqdotsc = [skewwsc * Mv;
              skewwsc * Mw;
   -(Lr'/m) * skewwsc * Mv;
              zeros(NXRW,1)];
    % vectorized comp of Cqdotrw: faster but less intuitive
    Cqdotrw = zeros(NXd,1);
    wrwi     = wsc' + tensorprod(Nrw, etadot,3,1);
    skewwrwi = skewVec(wrwi);
    Cqdotrw(IdXdrigr) = tensorprod(skewwrwi,Mwrw,[1 3],[1 2]);

    a = reshape(permute(skewwrwi,[2,1,3]),[],Nrigrd)*Mwrw';
    A = a(:);
    ijxi2ji=(0:NXRW-1)*((NXRW+1)*Nrigrd)+(1:Nrigrd)';
    ijxi2ji = ijxi2ji(:);
    b = A(ijxi2ji);
    Cqdotrw(IdXdmod) = reshape(permute(Nrw,[2,1,3]),[],Nmod)'*b;
    
    MinvCqdot = Minv * (Cqdotsc + Cqdotrw);

    if ~SimParam.SimRWPertFlag
        w     = zeros(IdStruct.NW,1);
        MinvG = zeros(NXd,IdStruct.NW);
    else
        % compute perturbations 
        [wKnown,wUnknown] = RWPertFunc(thetarw,wrw,ModelStruct,IdStruct,SimParam); 
        w = (wKnown + wUnknown)';
    end

    % filters
    % call meas function and apply filter
    if FiltFlag && ~isempty(IdXFilt)
        AFilt = ModelStruct.Sens.SysSensA;
        BFilt = ModelStruct.Sens.SysSensB;
        IdStruct.RealMdl.NFilt = 0;
        [~,~,~,meas_noerr] = ObsFunc(RealStruct, x, u, w, thetarw, ...
                                            IdStruct, ModelStruct, SimParam);
        
        xdot(IdXFilt) = AFilt * x(IdXFilt) + BFilt * meas_noerr;
    end
    
    % MRP prop
    xdot(IdStruct.RealMdl.IdXrigrq) = (pp/4*eye(3) + sp*(sp+eye(3))/2)*wsc;
    xdot(IdStruct.RealMdl.IdXmodq) = x(IdStruct.RealMdl.IdXmodd);
    xdot(IdStruct.RealMdl.IdXd)    = MinvB * u  + MinvG * w - MinvK * eta ...
                                             - MinvC * etadot - MinvCqdot;
    xdot(NX+1:end) = wrw;

end

