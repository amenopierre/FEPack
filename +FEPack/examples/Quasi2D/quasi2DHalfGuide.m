function [R, D, newBCstruct, Lambda] = quasi2DHalfGuide(meshXY, meshYZ, meshLineZ, orientation, FEcos, cutvec, BCstruct, quadratureRatio, suffix, folder_name, opts)

  %% Initialization
  cutslope = cutvec(2) / cutvec(1);
  N_s = meshLineZ.numPoints;
  % sVec = meshLineZ.points;
  N_XY = meshXY.numPoints;
  % cellXY = meshXY.domain('volumic');
  % cutvec = cutvec;
  % suffix = opts.suffix;
  % omega = opts.omega;
  fid = fopen('outputs/elapsed_time', 'a+');
  co_gradu_gradv = FEcos.gradu_gradv;
  co_gradu_vectv = FEcos.gradu_vectv;
  co_vectu_gradv = FEcos.vectu_gradv;
  co_vectu_vectv = FEcos.vectu_vectv;
  co_funQ_u_v    = FEcos.funQ_u_v;   
  co_funR_u_v    = FEcos.funR_u_v;   
  u = FEPack.pdes.PDEObject;
  v = dual(u);

  % Edges
  edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
  edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
  edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
  edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints; 
  
  % Surfacic rhs
  B0x = sparse(edge0x.IdPoints(1:end), (1:N0x), 1.0, N_XY, N0x);
  B1x = sparse(edge1x.IdPoints(1:end), (1:N1x), 1.0, N_XY, N1x);
  B0y = sparse(edge0y.IdPoints(2:end-1), (1:N0y-2), 1.0, N_XY, N0y-2);
  B1y = sparse(edge1y.IdPoints(2:end-1), (1:N1y-2), 1.0, N_XY, N1y-2);

  % Basis functions on face X = cst and S/Z-line  
  spBX = BCstruct.spBX; NbX = spBX.numBasis; % X = cst
  spBS = BCstruct.spBS; NbS = spBS.numBasis; % S/Z-line
  NbY = NbS * (N0y - 2);

  % Define shear maps
  % /////////////////
  % X = cst
  cut1_0x = cutvec(1) * meshXY.points(edge0x.IdPoints, 2);
  cut2_0x = cutvec(2) * meshXY.points(edge0x.IdPoints, 2);
  cut1_1x = cutvec(1) * meshXY.points(edge1x.IdPoints, 2);
  cut2_1x = cutvec(2) * meshXY.points(edge1x.IdPoints, 2);
  %
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([cut1_0x, FEPack.tools.mymod(cut2_0x + s, 0, 1), zeros(N0x, 1)]);
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([cut1_1x, FEPack.tools.mymod(cut2_1x + s, 0, 1), zeros(N1x, 1)]);
  %
  % S/Z-line
  [IdSX_S, IdSX_X] = ind2sub([NbS, N0y-2], 1:NbY);
  phisY = eye(N0y-2);
  shearmap.fun0y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPack.tools.mymod(s,            0, 1) 0, 0], IdSX_S);
  shearmap.fun1y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPack.tools.mymod(s + cutslope, 0, 1) 0, 0], IdSX_S);

  % Quadrature for DtN/RtR operators: trapezoidal rule
  Nquad = N_s * quadratureRatio;
  quadpoints = linspace(0, 1, Nquad + 1).';
  quadpoints(end) = [];
  structLoc = meshLineZ.domain('volumic').locateInDomain([quadpoints, zeros(Nquad, 2)]);
  elts = meshLineZ.domain('volumic').elements(structLoc.elements, :);
  coos = structLoc.barycoos;
  quadwgt = (quadpoints(2) - quadpoints(1)) * cutvec(1);
  

  if (abs(BCstruct.BCdu) < eps)

    % Dirichlet boundary conditions %
    % ***************************** %
    if (opts.verbose)
      fprintf('\tConditions de Dirichlet\n');
    end

    if (isfield(BCstruct, 'BCu') && ...
      (~isa(BCstruct.BCu, 'double') || (BCstruct.BCu ~= 1.0)))
      warning('on');
      warning(['Pour des conditions de Dirichlet, le coefficient devant u ',...
               'dans la condition aux bords est automatiquement pris égal à 1.']);
    end
    BCstruct.BCu = 1.0;
    BCstruct.representation = 'projection';

    % Homogeneous Dirichlet boundary conditions
    ecs = ((u|edge0x) == 0.0) & ((u|edge1x) == 0.0) &...
          ((u|edge0y) == 0.0) & ((u|edge1y) == 0.0);

    ecs.applyEcs;
    ecs.b = [B0x, B1x, B0y, B1y];
    sidenames = {'0x', '1x', '0y', '1y'};

    %% Local cell problems
    %  ///////////////////
    fprintf('<strong>Problèmes de cellule locaux.</strong>\n');
    tic;
    for idS = 1:N_s-1 % parfor
      
      % The FE matrix is a linear combination of elementary matrices
      FEmat = load([folder_name, 'FEmat_', suffix, '_', num2str(idS), '.mat']);

      AA = co_gradu_gradv * FEmat.mat_gradu_gradv...
         + co_gradu_vectv * FEmat.mat_gradu_vectv...
         + co_vectu_gradv * FEmat.mat_vectu_gradv...
         + co_vectu_vectv * FEmat.mat_vectu_vectv...
         + co_funQ_u_v    * FEmat.mat_funQ_u_v...
         + co_funR_u_v    * FEmat.mat_funR_u_v;
      
      AA0 =  ecs.P * AA * ecs.P'; %#ok
      LL0 = -ecs.P * AA * ecs.b;
      
      % Solve the linear system
      % tic;
      Ecell0 = ecs.b + ecs.P' * (AA0 \ LL0);
      % tpscpu = toc;
      % fprintf('Inversion système : %0.3e secondes\n', tpscpu);

      % Local cell solutions
      solcell = [];
      solcell.E0x = Ecell0(:, 1:N0x);   Ecell0(:, 1:N0x)   = [];
      solcell.E1x = Ecell0(:, 1:N1x);   Ecell0(:, 1:N1x)   = [];
      solcell.E0y = Ecell0(:, 1:N0y-2); Ecell0(:, 1:N0y-2) = [];
      solcell.E1y = Ecell0(:, 1:N1y-2);

      % DtN operators
      for idI = 1:4
        nameEsolI = ['E', sidenames{idI}]; %#ok
        AA_m_solI = AA * solcell.(nameEsolI);

        for idJ = 1:4
          nameEsolJ = ['E', sidenames{idJ}];
          nameEdgeT = ['edgeT', sidenames{idI}, sidenames{idJ}];
          solcell.(nameEdgeT) = solcell.(nameEsolJ)' * AA_m_solI;
        end
      end

      % Save local cell solutions
      FEPack.tools.parsave([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(idS)], solcell, true);
    end
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);
    system(['cp ', folder_name, 'local_cell_sol_', suffix, '_S_1.mat ', folder_name, 'local_cell_sol_', suffix, '_S_', num2str(N_s), '.mat']);
    
    %% DtN operators
    %  /////////////
    % Initialize the Face DtN operators
    fprintf('<strong>Opérateurs DtN.</strong>\n');
    T0x0x = zeros(NbX, NbX); T0x1x = T0x0x; T1x0x = T0x0x; T1x1x = T0x0x; 
    T0x0y = zeros(NbY, NbX); T0x1y = T0x0y; T1x0y = T0x0y; T1x1y = T0x0y; 
    T0y0x = zeros(NbX, NbY); T0y1x = T0y0x; T1y0x = T0y0x; T1y1x = T0y0x; 
    T0y0y = zeros(NbY, NbY); T0y1y = T0y0y; T1y0y = T0y0y; T1y1y = T0y0y;
    
    % Local auxiliary DtN operators
    tic;
    for idQuad = 1:Nquad % parfor

      % Local auxiliary DtN operators
      sVar = quadpoints(idQuad);
      eltsQuad = elts(idQuad, :);
      coosQuad = coos(idQuad, :);
      sc1 = load([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(eltsQuad(1))]);
      sc2 = load([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(eltsQuad(2))]);

      T0x0x = T0x0x + shearmap.fun0x(sVar)' * (coosQuad(1) * sc1.edgeT0x0x + coosQuad(2) * sc2.edgeT0x0x) * shearmap.fun0x(sVar); %#ok 
      T1x0x = T1x0x + shearmap.fun0x(sVar)' * (coosQuad(1) * sc1.edgeT1x0x + coosQuad(2) * sc2.edgeT1x0x) * shearmap.fun1x(sVar); 
      T0x1x = T0x1x + shearmap.fun1x(sVar)' * (coosQuad(1) * sc1.edgeT0x1x + coosQuad(2) * sc2.edgeT0x1x) * shearmap.fun0x(sVar); 
      T1x1x = T1x1x + shearmap.fun1x(sVar)' * (coosQuad(1) * sc1.edgeT1x1x + coosQuad(2) * sc2.edgeT1x1x) * shearmap.fun1x(sVar);
      
      T0x0y = T0x0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT0x0y + coosQuad(2) * sc2.edgeT0x0y) * shearmap.fun0x(sVar); 
      T1x0y = T1x0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT1x0y + coosQuad(2) * sc2.edgeT1x0y) * shearmap.fun1x(sVar); 
      T0x1y = T0x1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT0x1y + coosQuad(2) * sc2.edgeT0x1y) * shearmap.fun0x(sVar); 
      T1x1y = T1x1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT1x1y + coosQuad(2) * sc2.edgeT1x1y) * shearmap.fun1x(sVar);
      
      T0y0x = T0y0x + shearmap.fun0x(sVar)' * (coosQuad(1) * sc1.edgeT0y0x + coosQuad(2) * sc2.edgeT0y0x) * shearmap.fun0y(sVar); 
      T1y0x = T1y0x + shearmap.fun0x(sVar)' * (coosQuad(1) * sc1.edgeT1y0x + coosQuad(2) * sc2.edgeT1y0x) * shearmap.fun1y(sVar); 
      T0y1x = T0y1x + shearmap.fun1x(sVar)' * (coosQuad(1) * sc1.edgeT0y1x + coosQuad(2) * sc2.edgeT0y1x) * shearmap.fun0y(sVar); 
      T1y1x = T1y1x + shearmap.fun1x(sVar)' * (coosQuad(1) * sc1.edgeT1y1x + coosQuad(2) * sc2.edgeT1y1x) * shearmap.fun1y(sVar);
      
      T0y0y = T0y0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT0y0y + coosQuad(2) * sc2.edgeT0y0y) * shearmap.fun0y(sVar); 
      T1y0y = T1y0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT1y0y + coosQuad(2) * sc2.edgeT1y0y) * shearmap.fun1y(sVar); 
      T0y1y = T0y1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT0y1y + coosQuad(2) * sc2.edgeT0y1y) * shearmap.fun0y(sVar); 
      T1y1y = T1y1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT1y1y + coosQuad(2) * sc2.edgeT1y1y) * shearmap.fun1y(sVar);

    end
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);

    %% Compute DtD operator
    %  ////////////////////
    fprintf('<strong>Calcul opérateur DtD.</strong>\n');

    % Divide DtN operators by the quadrature weight
    T0x0x = T0x0x * quadwgt;
    T1x0x = T1x0x * quadwgt;
    T0x1x = T0x1x * quadwgt;
    T1x1x = T1x1x * quadwgt;
    T0x0y = T0x0y * quadwgt;
    T1x0y = T1x0y * quadwgt;
    T0x1y = T0x1y * quadwgt;
    T1x1y = T1x1y * quadwgt;
    T0y0x = T0y0x * quadwgt;
    T1y0x = T1y0x * quadwgt;
    T0y1x = T0y1x * quadwgt;
    T1y1x = T1y1x * quadwgt;
    T0y0y = T0y0y * quadwgt;
    T1y0y = T1y0y * quadwgt;
    T0y1y = T0y1y * quadwgt;
    T1y1y = T1y1y * quadwgt;

    % Solve linear system
    tic;
    DtD = (T1y0y + T0y0y + T1y1y + T0y1y) \ [-(T0x0y + T0x1y), -(T1x0y + T1x1y)];
    solguide.DtD0 = DtD(:,     1:NbX);
    solguide.DtD1 = DtD(:, NbX+1:end);
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);

    %% Traces and normal traces of the local cell solutions
    %  ////////////////////////////////////////////////////
    fprintf('<strong>Opérateurs DtN locaux.</strong>\n');

    % Ekl is the trace of Ek on Sigmal
    E00 = speye(NbX);
    E10 = sparse(NbX, NbX);
    E01 = sparse(NbX, NbX);
    E11 = speye(NbX);

    % Fkl is the normal trace of Ek on Sigmal. It is evaluated weakly and then projected
    F00 = spBX.invmassmat * ((T0y0x + T1y0x) * solguide.DtD0 + T0x0x);
    F10 = spBX.invmassmat * ((T0y0x + T1y0x) * solguide.DtD1 + T1x0x);
    F01 = spBX.invmassmat * ((T0y1x + T1y1x) * solguide.DtD0 + T0x1x);
    F11 = spBX.invmassmat * ((T0y1x + T1y1x) * solguide.DtD1 + T1x1x);

  else

    % Robin (Neumann) boundary conditions %
    % *********************************** %
    if (opts.verbose)
      fprintf('\tConditions de Robin-Neumann\n');
    end

    % Homogeneous Dirichlet boundary conditions on y = cst only
    ecs = ((u|edge0y) == 0.0) & ((u|edge1y) == 0.0);

    ecs.applyEcs;
    ecs.b = [zeros(N_XY, N0x + N1x), B0y, B1y];
    sidenames = {'0x', '1x', '0y', '1y'};

    % Matrix associated to the multiplication by a function
    if isa(BCstruct.BCu, 'function_handle')
      BCstruct.BCu = spB0.intg_U_V(edge0x, BCstruct.BCu);
      BCstruct.representation = 'weak evaluation';
    end

    % Surfacic contributions
    if (~isfield(BCstruct, 'representation') || isscalar(BCstruct.BCu))
      BCstruct.representation = 'projection';
    end

    SS0 = FEPack.pdes.Form.intg_TU_V(edge0x, BCstruct.BCu, BCstruct.representation);
    SS1 = FEPack.pdes.Form.intg_TU_V(edge1x, BCstruct.BCu, BCstruct.representation);

    % Surfacic right-hand sides
    LL = [FEPack.pdes.Form.intg(edge0x, u*v) * B0x,...
          FEPack.pdes.Form.intg(edge1x, u*v) * B1x, zeros(N_XY, N0y+N1y-4)];
    
    % Auxiliary variables
    invBCdu = 1.0 / BCstruct.BCdu;
    BCu = BCstruct.BCu;
    edge0x_IdPoints = edge0x.IdPoints;
    edge1x_IdPoints = edge1x.IdPoints;
    speye_y0x = speye(N0x); speye_y0x(:, 1) = []; speye_y0x(:, end) = [];
    speye_y1x = speye(N1x); speye_y1x(:, 1) = []; speye_y1x(:, end) = [];

    %% Local cell problems
    %  ///////////////////
    fprintf('<strong>Problèmes de cellule locaux.</strong>\n');
    tic;
    for idS = 1:N_s-1 % parfor
      
      % The FE matrix is a linear combination of elementary matrices
      FEmat = load([folder_name, 'FEmat_', suffix, '_', num2str(idS), '.mat']);

      AA = co_gradu_gradv * FEmat.mat_gradu_gradv...
         + co_gradu_vectv * FEmat.mat_gradu_vectv...
         + co_vectu_gradv * FEmat.mat_vectu_gradv...
         + co_vectu_vectv * FEmat.mat_vectu_vectv...
         + co_funQ_u_v    * FEmat.mat_funQ_u_v...
         + co_funR_u_v    * FEmat.mat_funR_u_v;
      
      AA = AA + SS0 + SS1;
      
      AA0 =  ecs.P * AA * ecs.P'; %#ok
      LL0 = -ecs.P * AA * ecs.b + ecs.P * LL;

      % Solve the linear system
      % tic;
      Ecell0 = ecs.b + ecs.P' * (AA0 \ LL0);
      % tpscpu = toc;
      % fprintf('Inversion système : %0.3e secondes\n', tpscpu);

      % Local cell solutions
      solcell = [];
      solcell.E0x = Ecell0(:, 1:N0x);   Ecell0(:, 1:N0x)   = [];
      solcell.E1x = Ecell0(:, 1:N1x);   Ecell0(:, 1:N1x)   = [];
      solcell.E0y = Ecell0(:, 1:N0y-2); Ecell0(:, 1:N0y-2) = [];
      solcell.E1y = Ecell0(:, 1:N1y-2);

      % Edge transmission operators
      % === Operators associated to Neumann trace on y = cst only
      for idI = 1:4
        nameEsolI = ['E', sidenames{idI}]; %#ok
        AA_m_solI = AA * solcell.(nameEsolI);

        for idJ = 3:4 % For the trace on y = cst
          nameEsolJ = ['E', sidenames{idJ}];
          nameEdgeT = ['edgeT', sidenames{idI}, sidenames{idJ}];
          solcell.(nameEdgeT) = solcell.(nameEsolJ)' * AA_m_solI;
        end
      end

      % === Operators associated to Robin trace on x = cst only
      solcell.edgeT0x0x = -invBCdu * (BCu * solcell.E0x(edge0x_IdPoints, :) - speye(N0x));
      solcell.edgeT1x0x = -invBCdu *  BCu * solcell.E1x(edge0x_IdPoints, :);
      solcell.edgeT0x1x = -invBCdu *  BCu * solcell.E0x(edge1x_IdPoints, :);
      solcell.edgeT1x1x = -invBCdu * (BCu * solcell.E1x(edge1x_IdPoints, :) - speye(N1x));

      solcell.edgeT0y0x = -invBCdu * (BCu * solcell.E0y(edge0x_IdPoints, :) - speye_y0x);
      solcell.edgeT1y0x = -invBCdu *  BCu * solcell.E1y(edge0x_IdPoints, :);
      solcell.edgeT0y1x = -invBCdu *  BCu * solcell.E0y(edge1x_IdPoints, :);
      solcell.edgeT1y1x = -invBCdu * (BCu * solcell.E1y(edge1x_IdPoints, :) - speye_y1x);

      % Save local cell solutions
      FEPack.tools.parsave([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(idS)], solcell, true);
    end
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);
    system(['cp ', folder_name, 'local_cell_sol_', suffix, '_S_1.mat ', folder_name, 'local_cell_sol_', suffix, '_S_', num2str(N_s), '.mat']);

    %% Transmission operators
    % ///////////////////////
    % Initialize the Face DtN operators
    fprintf('<strong>Opérateurs DtN.</strong>\n');
    T0x0y = zeros(NbY, NbX); T0x1y = T0x0y; T1x0y = T0x0y; T1x1y = T0x0y; 
    T0y0y = zeros(NbY, NbY); T0y1y = T0y0y; T1y0y = T0y0y; T1y1y = T0y0y;
    
    % Local auxiliary DtN operators for y = cst only
    tic;
    for idQuad = 1:Nquad % parfor

      % Local auxiliary DtN operators
      sVar = quadpoints(idQuad);
      eltsQuad = elts(idQuad, :);
      coosQuad = coos(idQuad, :);
      sc1 = load([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(eltsQuad(1))]);
      sc2 = load([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(eltsQuad(2))]);

      T0x0y = T0x0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT0x0y + coosQuad(2) * sc2.edgeT0x0y) * shearmap.fun0x(sVar); %#ok
      T1x0y = T1x0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT1x0y + coosQuad(2) * sc2.edgeT1x0y) * shearmap.fun1x(sVar); 
      T0x1y = T0x1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT0x1y + coosQuad(2) * sc2.edgeT0x1y) * shearmap.fun0x(sVar); 
      T1x1y = T1x1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT1x1y + coosQuad(2) * sc2.edgeT1x1y) * shearmap.fun1x(sVar);
      
      T0y0y = T0y0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT0y0y + coosQuad(2) * sc2.edgeT0y0y) * shearmap.fun0y(sVar); 
      T1y0y = T1y0y + shearmap.fun0y(sVar)' * (coosQuad(1) * sc1.edgeT1y0y + coosQuad(2) * sc2.edgeT1y0y) * shearmap.fun1y(sVar); 
      T0y1y = T0y1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT0y1y + coosQuad(2) * sc2.edgeT0y1y) * shearmap.fun0y(sVar); 
      T1y1y = T1y1y + shearmap.fun1y(sVar)' * (coosQuad(1) * sc1.edgeT1y1y + coosQuad(2) * sc2.edgeT1y1y) * shearmap.fun1y(sVar);

    end
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);
    
    % Local auxiliary RtR operators for x = cst only
    fprintf('<strong>Opérateurs RtR.</strong>\n');
    T0x0x = zeros(NbX, NbX); T0x1x = T0x0x; T1x0x = T0x0x; T1x1x = T0x0x; 
    T0y0x = zeros(NbX, NbY); T0y1x = T0y0x; T1y0x = T0y0x; T1y1x = T0y0x; 
    
    spts = FEPack.tools.mymod(meshYZ.points(:, 2) - meshYZ.points(:, 1) * cutslope, 0, 1);
    Xpts = meshYZ.points(:, 1) / cutvec(1); 
    N_YZ = meshYZ.numPoints;

    structLoc_s = meshLineZ.domain('volumic').locateInDomain([spts, zeros(N_YZ, 2)]);
    elts_spts = meshLineZ.domain('volumic').elements(structLoc_s.elements, :);
    coos_spts = structLoc_s.barycoos;

    structLoc0x = edge0x.locateInDomain([Xpts, zeros(N_YZ, 2)]);
    elts0x = meshLineZ.domain('volumic').elements(structLoc0x.elements, :);
    coos0x = structLoc0x.barycoos;

    structLoc1x = edge1x.locateInDomain([Xpts, zeros(N_YZ, 2)]);
    elts1x = meshLineZ.domain('volumic').elements(structLoc1x.elements, :);
    coos1x = structLoc1x.barycoos;
    
    tic;
    for idYZ = 1:N_YZ % parfor

      % Local auxiliary DtN operators
      sVar = spts(idYZ);
      elts_s = elts_spts(idYZ, :);
      coos_s = coos_spts(idYZ, :);
      sc1 = load([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(elts_s(1))]);
      sc2 = load([folder_name, 'local_cell_sol_', suffix, '_S_', num2str(elts_s(2))]);

      elt0x = elts0x(idYZ, :); elt1x = elts1x(idYZ, :);
      coo0x = coos0x(idYZ, :); coo1x = coos1x(idYZ, :);

      t0x0x_on_edge = (coos_s(1) * sc1.edgeT0x0x + coos_s(2) * sc2.edgeT0x0x) * shearmap.fun0x(sVar); %#ok
      t1x0x_on_edge = (coos_s(1) * sc1.edgeT1x0x + coos_s(2) * sc2.edgeT1x0x) * shearmap.fun0x(sVar);
      t0x1x_on_edge = (coos_s(1) * sc1.edgeT0x1x + coos_s(2) * sc2.edgeT0x1x) * shearmap.fun1x(sVar);
      t1x1x_on_edge = (coos_s(1) * sc1.edgeT1x1x + coos_s(2) * sc2.edgeT1x1x) * shearmap.fun1x(sVar);

      t0y0x_on_edge = (coos_s(1) * sc1.edgeT0y0x + coos_s(2) * sc2.edgeT0y0x) * shearmap.fun0y(sVar);
      t1y0x_on_edge = (coos_s(1) * sc1.edgeT1y0x + coos_s(2) * sc2.edgeT1y0x) * shearmap.fun0y(sVar);
      t0y1x_on_edge = (coos_s(1) * sc1.edgeT0y1x + coos_s(2) * sc2.edgeT0y1x) * shearmap.fun1y(sVar);
      t1y1x_on_edge = (coos_s(1) * sc1.edgeT1y1x + coos_s(2) * sc2.edgeT1y1x) * shearmap.fun1y(sVar);


      T0x0x(idYZ, :) = coo0x(1) * t0x0x_on_edge(elt0x(1), :) + coo0x(2) * t0x0x_on_edge(elt0x(2), :);
      T1x0x(idYZ, :) = coo0x(1) * t1x0x_on_edge(elt0x(1), :) + coo0x(2) * t1x0x_on_edge(elt0x(2), :);
      T0x1x(idYZ, :) = coo1x(1) * t0x1x_on_edge(elt1x(1), :) + coo1x(2) * t0x1x_on_edge(elt1x(2), :);
      T1x1x(idYZ, :) = coo1x(1) * t1x1x_on_edge(elt1x(1), :) + coo1x(2) * t1x1x_on_edge(elt1x(2), :);
      
      T0y0x(idYZ, :) = coo0x(1) * t0y0x_on_edge(elt0x(1), :) + coo0x(2) * t0y0x_on_edge(elt0x(2), :);
      T1y0x(idYZ, :) = coo0x(1) * t1y0x_on_edge(elt0x(1), :) + coo0x(2) * t1y0x_on_edge(elt0x(2), :);
      T0y1x(idYZ, :) = coo1x(1) * t0y1x_on_edge(elt1x(1), :) + coo1x(2) * t0y1x_on_edge(elt1x(2), :);
      T1y1x(idYZ, :) = coo1x(1) * t1y1x_on_edge(elt1x(1), :) + coo1x(2) * t1y1x_on_edge(elt1x(2), :);

    end
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);

    %% Compute DtD operator
    %  ////////////////////
    fprintf('<strong>Calcul opérateur DtD.</strong>\n');

    % Divide DtN operators by the quadrature weight
    T0x0y = T0x0y * quadwgt;
    T1x0y = T1x0y * quadwgt;
    T0x1y = T0x1y * quadwgt;
    T1x1y = T1x1y * quadwgt;
    T0y0y = T0y0y * quadwgt;
    T1y0y = T1y0y * quadwgt;
    T0y1y = T0y1y * quadwgt;
    T1y1y = T1y1y * quadwgt;

    % Solve linear system
    tic;
    DtD = (T1y0y + T0y0y + T1y1y + T0y1y) \ [-(T0x0y + T0x1y), -(T1x0y + T1x1y)];
    solguide.DtD0 = DtD(:,     1:NbX);
    solguide.DtD1 = DtD(:, NbX+1:end);
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);

    %% Traces and normal traces of the local cell solutions
    %  ////////////////////////////////////////////////////
    fprintf('<strong>Opérateurs DtN locaux.</strong>\n');

    % Ekl is the trace of Ek on Sigmal
    E00 = speye(NbX);
    E10 = sparse(NbX, NbX);
    E01 = sparse(NbX, NbX);
    E11 = speye(NbX);

    % Fkl is the normal trace of Ek on Sigmal. It is evaluated weakly and then projected
    F00 = spBX.projmat * ((T0y0x + T1y0x) * solguide.DtD0 + T0x0x);
    F10 = spBX.projmat * ((T0y0x + T1y0x) * solguide.DtD1 + T1x0x);
    F01 = spBX.projmat * ((T0y1x + T1y1x) * solguide.DtD0 + T0x1x);
    F11 = spBX.projmat * ((T0y1x + T1y1x) * solguide.DtD1 + T1x1x);

  end

  %% Solve the Riccati equation
  fprintf('<strong>Equation de Riccati.</strong>\n');
  
  flux = @(V) modesFlux(V, E00, E10, F00, F10, orientation, spBX, opts.omega);
  riccatiOpts.tol = 1.0e-2;
  if isfield(opts, 'suffix')
    riccatiOpts.suffix = suffix;
  else
    riccatiOpts.suffix = '';
  end
  % riccatiOpts.suffix = [riccatiOpts.suffix, '_', int2str(idFB)];

  % Solve the linearized eigenvalue problem associated to the Riccati equation
  tic;
  [R, D] = propagationOperators([E01, E11;  orientation*F01,  orientation*F11], ...
                                [E00, E10; -orientation*F00, -orientation*F10], flux, riccatiOpts);
  tps = toc;
  fprintf(fid, '%0.5e\n', tps);
  % fprintf('Résolution équation de Riccati : %0.3e secondes\n', tps);
  
  % % Cas homogène
  % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FourierIdsX = spBX.FourierIds.X(:).'; dx = length(FourierIdsX);
  % FourierIdsY = spBX.FourierIds.Y(:).'; dy = length(FourierIdsY);
  % FourierIdsZ = spBX.FourierIds.Z(:).'; dz = length(FourierIdsZ);

  % [Ix, Iy, ~] = ind2sub([dx, dy, dz], (1:spBX.numBasis)');
  % r_xi_fun = FloquetVar*cutvec(1) + 2*pi*(FourierIdsX(Ix) * cutvec(1) + FourierIdsY(Iy) * cutvec(2)).';
  % l_xi_fun = sqrt(r_xi_fun.^2 - opts.omega^2);

  % % mZXpts = spBY.domain.mesh.points;

  % % solguide.DtD0 = spBY.FE_to_spectral * (sinh((1 - orientation*mZXpts(:, 2)) * l_xi_fun.') .* exp(2i*pi*mZXpts(:, 1) * FourierIdsY(Iy)) ./ (ones(size(mZXpts, 1), 1) * sinh(l_xi_fun.')));
  % % %
  % % % norm(solguide.DtD0)
  % % %
  % % solguide.DtD1 = spBY.FE_to_spectral * (sinh(orientation*mZXpts(:, 2) * l_xi_fun.') .* exp(2i*pi*mZXpts(:, 1) * FourierIdsY(Iy)) ./ (ones(size(mZXpts, 1), 1) * sinh(l_xi_fun.')));

  % F00 =  diag(l_xi_fun./tanh(l_xi_fun));
  % F01 = -diag(l_xi_fun./sinh(l_xi_fun));
  % F10 = F01;
  % F11 = F00;

  % solguide.R = diag(exp(-l_xi_fun));
  % solguide.D = solguide.R;

  % E00 = speye(NbX);
  % E10 = sparse(NbX, NbX);
  % E01 = sparse(NbX, NbX);
  % E11 = speye(NbX);
  % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% Compute the DtN/RtR operator
  U0 = E00 + E10 * D;
  dU0 = F00 + F10 * D;

  Lambda = -BCstruct.BCu' * dU0 + BCstruct.BCdu' * U0;
  newBCstruct = BCstruct;

  fclose(fid);

end

%% Auxiliary function: The flux function
function flux = modesFlux(V, E00, E10, F00, F10, orientation, spB, omega)
  
  Nb = size(E00, 1);
  Psi  = E00 * V(1:Nb, :) + E10 * V(Nb+1:end, :);
  dPsi = F00 * V(1:Nb, :) + F10 * V(Nb+1:end, :);

  flux = -orientation * imag(diag(Psi' * spB.massmat * dPsi) / omega);

end
