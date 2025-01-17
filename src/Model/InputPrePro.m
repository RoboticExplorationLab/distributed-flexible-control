function [ModelStruct,IdStruct] = InputPrePro(TableFile, SimParam)
    ModelStruct = struct;

    Table = TableFile.T;
    % FEM Mass Matrix
    MassFEM = Table.M{:};
    % FEM Stiffness Matrix
    StiffFEM = Table.K{:};
    % Connectivity Matrix
    ConnectMat = Table.ConnectMat{:};
    % coordinates of the nodes in the FEM mesh
    NodeCoord = Table.NodeCoord{:};
    % Indices of the boundary nodes
    IdBound = Table.IdBound{:};
    % structure with FEM data
    
    %% Input pre-processing
    % FEM Mapping Matrix: node to element coordinates.
    Nfuncs = Table.N{:};
    
    %% Input pre-processing
    % mass matrix of the antenna needs to be shifted according to its position
    % wrt the SC
    % boundary nodes
    MassFEM(9+IdBound,:) = [];
    MassFEM(:,9+IdBound) = [];
    % additional coords irrelevant to this project
    MassFEM(7:9,:)       = [];
    MassFEM(:,7:9)       = [];
    
    % Stiffness matrix K
    % StiffFEM = Table.K{:}/(2^4); % /(2^2); % ; % /(2^5); % 4; % 
    StiffFEM = Table.K{:}/(2^6); % /(2^2); % ; % /(2^5); % 4; % 
    StiffFEM(IdBound,:) = [];
    StiffFEM(:,IdBound) = [];

    Nnode     = size(NodeCoord,1);
    IdBound = unique(ceil(IdBound/3));
    IdInner   = setdiff(1:Nnode,IdBound);

    ModelStruct.Nfuncs = Nfuncs;

    %% State indexing
    % Structure of indexes
    IdStruct           = struct;
    % rigid translational states are propagated in nonlinear model, but
    % removed from the linear one

    % Structure of indexes for true nonlinear model
    IdStruct.RealMdl = struct;

    % number and Id of rigid body pos coord
    IdStruct.RealMdl.Nrigq  = 6;
    % number and Id of rigid body pos trans coord
    IdStruct.RealMdl.Nrigtq = 3;
    % number and Id of rigid body pos rot coord
    IdStruct.RealMdl.Nrigrq = 3;
    % number and Id of rigid body vel coord
    IdStruct.RealMdl.Nrigd  = 6;
    % number and Id of rigid body vel trans coord
    IdStruct.RealMdl.Nrigtd = 3;
    % number and Id of rigid body vel rot coord
    IdStruct.RealMdl.Nrigrd = 3;

    % number and Id of FEM coordinates
    IdStruct.RealMdl.Nfem  = size(MassFEM,1) - IdStruct.RealMdl.Nrigd; 

    % total number and Id of rigid body coordinates
    IdStruct.RealMdl.Nrig  = IdStruct.RealMdl.Nrigq ...
                                          + IdStruct.RealMdl.Nrigd;
    IdStruct.RealMdl.Nrigt = IdStruct.RealMdl.Nrigtq ...
                                         + IdStruct.RealMdl.Nrigtd;
    IdStruct.RealMdl.Nrigr = IdStruct.RealMdl.Nrigrq ...
                                         + IdStruct.RealMdl.Nrigrd;
    
    % % total number of position states and indices within them
    IdStruct.RealMdl.NXq      = IdStruct.RealMdl.Nrigq;
    IdStruct.RealMdl.IdXq     = 1:IdStruct.RealMdl.NXq;
    % index within pos states
    IdStruct.RealMdl.IdXqrig  = 1:IdStruct.RealMdl.Nrigq;
    IdStruct.RealMdl.IdXqrigt = 1:IdStruct.RealMdl.Nrigtq;
    IdStruct.RealMdl.IdXqrigr = IdStruct.RealMdl.Nrigtq ...
                                   + (1:IdStruct.RealMdl.Nrigrq);
    % index within all states
    IdStruct.RealMdl.IdXrigq  = IdStruct.RealMdl.IdXqrig;
    IdStruct.RealMdl.IdXrigtq = IdStruct.RealMdl.IdXqrigt;
    IdStruct.RealMdl.IdXrigrq = IdStruct.RealMdl.IdXqrigr;
    
    % % total number of velocity states and indices within them
    IdStruct.RealMdl.NXd      = IdStruct.RealMdl.Nrigd;
    IdStruct.RealMdl.IdXd     = IdStruct.RealMdl.NXq ...
                                        + (1:IdStruct.RealMdl.NXd);
    % index within vel states
    IdStruct.RealMdl.IdXdrig  = 1:IdStruct.RealMdl.Nrigd;
    IdStruct.RealMdl.IdXdrigt = 1:IdStruct.RealMdl.Nrigtd;
    IdStruct.RealMdl.IdXdrigr = IdStruct.RealMdl.Nrigtd ...
                                    + (1:IdStruct.RealMdl.Nrigrd);
    % index within all states
    IdStruct.RealMdl.IdXrigd  = IdStruct.RealMdl.Nrigq ...
                                       + IdStruct.RealMdl.IdXdrig;
    IdStruct.RealMdl.IdXrigtd = IdStruct.RealMdl.Nrigq ...
                                       + IdStruct.RealMdl.IdXdrigt;
    IdStruct.RealMdl.IdXrigrd = IdStruct.RealMdl.Nrigq ...
                                       + IdStruct.RealMdl.IdXdrigr;
    
    % total number of states
    IdStruct.RealMdl.NX      = IdStruct.RealMdl.NXq ...
                                            + IdStruct.RealMdl.NXd;
    IdStruct.RealMdl.IdXrig  = [IdStruct.RealMdl.IdXrigq, ...
                                        IdStruct.RealMdl.IdXrigd];

    % % save vars
    ModelStruct.ConnectMat = ConnectMat;
    ModelStruct.NodeCoord  = NodeCoord;
    ModelStruct.IdBound    = IdBound;
    ModelStruct.MassFEM    = MassFEM;
    ModelStruct.StiffFEM   = StiffFEM;
    ModelStruct.IdInner    = IdInner;
end

