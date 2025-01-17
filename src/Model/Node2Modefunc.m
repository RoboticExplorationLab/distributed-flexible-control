function [ModelStruct, IdStruct] = Node2Modefunc(ModelStruct, IdStruct, SimParam)
    % % load var
    MassFEM   = ModelStruct.MassFEM;
    StiffFEM  = ModelStruct.StiffFEM;
    RealId    = IdStruct.RealMdl;
    %% State indexing
    % number of vibration modes on the sim
    RealId.Nmod     = SimParam.Nmodes;

    % total number and Id of rigid body coordinates
    RealId.IdXmodq  = RealId.Nrigq + (1:RealId.Nmod);
    RealId.IdXmodd  = RealId.Nrig + RealId.Nmod + (1:RealId.Nmod);
    RealId.IdXmod   = [RealId.IdXmodq, RealId.IdXmodd];
    
    % total number of position states and indices within them
    RealId.NXq      = RealId.Nrigq + RealId.Nmod;
    RealId.IdXq     = 1:RealId.NXq;
    RealId.IdXqmod  = RealId.Nrigq + (1:RealId.Nmod);
    
    % total number of velocity states and indices within them
    RealId.NXd      = RealId.Nrigd + RealId.Nmod;
    RealId.IdXd     = RealId.NXq + (1:RealId.NXd);
    RealId.IdXrigd  = RealId.Nrigq + RealId.Nmod + RealId.IdXdrig;
    RealId.IdXrigtd = RealId.Nrigq + RealId.Nmod + RealId.IdXdrigt;
    RealId.IdXrigrd = RealId.Nrigq + RealId.Nmod + RealId.IdXdrigr;
    RealId.IdXdmod  = RealId.Nrigd + (1:RealId.Nmod);
    
    % total number of states
    RealId.NX      = RealId.NXq + RealId.NXd;
    RealId.IdXrig  = [RealId.IdXrigq, RealId.IdXrigd];
    
    %% Change of coordinates
    Mii = MassFEM(RealId.Nrigd+1:end, RealId.Nrigd+1:end); % fem only

    % decoupling rigid translational states
    FEMPartTransMat = MassFEM(1:RealId.Nrigtd, RealId.Nrigd+1:end);
    FEMPartRotMat   = MassFEM(RealId.IdXdrigr, RealId.Nrigd+1:end);

    [V,D] = eig(Mii\StiffFEM);
    [MinFreqs,IdFreqRob] = mink(diag(D),RealId.Nfem,1);
    
    % Mass-normalized Mode shapes 
    % q = Phi*eta, eta being the mode shape coordinates and q the nodal 
    % coordinates
    a         = chol(V'*Mii*V);
    Node2Mode = V/diag(diag(a));
    Node2Mode = Node2Mode(:,IdFreqRob);
    Mode2Node = Node2Mode\eye(size(Node2Mode,1));

    % Stiffness matrix for mode shape coordinates (diag(wj^2), wj nat. freq.)
    StiffMatModes = D(IdFreqRob,IdFreqRob);

    % Rotational modal participation matrix
    FEMPartRotMat   = MassFEM(RealId.IdXdrigr, RealId.Nrigd+1:end);
    ModPartRotMat   = FEMPartRotMat*Node2Mode;
    ModPartTransMat = FEMPartTransMat*Node2Mode;
    SCInertia       = MassFEM(RealId.IdXdrigr, RealId.IdXdrigr);
    SCMass          = MassFEM(1,1);
    % Psi = -K_b\Kij;

    % Modal Truncation (lowest modes)
    % SimParam.modeRng
    % find id of modes within this range
    freqs = (diag(StiffMatModes).^0.5)/(2*pi);
    Id    = find(and(freqs > SimParam.modeRng(1),freqs < SimParam.modeRng(2)));

    IdModes = Id(1:SimParam.Nmodes);
    ModPartRotMat   = ModPartRotMat(:,IdModes);
    ModPartTransMat = ModPartTransMat(:,IdModes);
    Node2Mode = Node2Mode(:,IdModes);
    Mode2Node = Mode2Node(IdModes,:);

    StiffMatModes = StiffMatModes(IdModes,IdModes);
    % Mass Matrix for mode shape coordinates (identity matrix)
    MassMatModes = eye(SimParam.Nmodes);
    
    % can be computed to keep damping ratios realistic within a given band
    % band of interest: lowest mode to nyquist frequency
    beta = SimParam.beta;
    DampMatModes = beta(1)*StiffMatModes + beta(2)*MassMatModes;
    % DampMatModes = 2*0.005*StiffMatModes.^0.5;
    
    % Descriptor matrix of full linear system
    % Exdot = Ax+Bu
    DescMatXq = eye(RealId.NXq);
    % Full Mass Matrix/ Descriptor Matrix for the velocity states
    DescMatXd = zeros(RealId.NXd);
    DescMatXd(RealId.IdXdrigt, RealId.IdXdrigt) = SCMass * eye(3);
    DescMatXd(RealId.IdXdrigr, RealId.IdXdrigr) = SCInertia;
    DescMatXd(RealId.IdXdrigt, RealId.IdXdmod)  = ModPartTransMat;
    DescMatXd(RealId.IdXdrigr, RealId.IdXdmod)  = ModPartRotMat;
    DescMatXd(RealId.IdXdmod, RealId.IdXdrigt)  = ModPartTransMat';
    DescMatXd(RealId.IdXdmod, RealId.IdXdrigr)  = ModPartRotMat';
    DescMatXd(RealId.IdXdmod, RealId.IdXdmod)   = MassMatModes;
    DescMat = blkdiag(DescMatXq,DescMatXd);  % sparse();
    
    % State Matrix A (from descriptor system)
    StateMat = zeros(RealId.NX);
    StateMat(RealId.IdXmodq, RealId.IdXmodd) = eye(SimParam.Nmodes);
    StateMat(RealId.IdXmodd, RealId.IdXmodq) = -StiffMatModes;
    StateMat(RealId.IdXmodd, RealId.IdXmodd) = -DampMatModes;
    
    % integration of rigid body motion
    StateMat(RealId.IdXrigtq, RealId.IdXrigtd) = eye(3);
    StateMat(RealId.IdXrigrq, RealId.IdXrigrd) = eye(3)/4;
    
    % % save vars
    ModelStruct.DescMat       = DescMat;
    ModelStruct.DampMatModes  = DampMatModes;
    ModelStruct.MassMatModes  = MassMatModes;
    ModelStruct.StiffMatModes = StiffMatModes;
    ModelStruct.Mode2Node     = Mode2Node;
    ModelStruct.Node2Mode     = Node2Mode;
    ModelStruct.StateMat      = StateMat;
    IdStruct.RealMdl          = RealId;
    ModelStruct.ModPartTransMat = ModPartTransMat;
    ModelStruct.ModPartRotMat   = ModPartRotMat;
    ModelStruct.SCInertia       = SCInertia;
    ModelStruct.SCMass          = SCMass;
end
