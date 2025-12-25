function quasi2DHalfGuideDirichlet_par(orientation, meshXY, meshLineZ, mu3D, rho3D, BCstruct, opts)

  %% Initialization
  cutslope = opts.cutvec(2) / opts.cutvec(1);
  nomdossier = opts.nomdossier;

  % Integrands for FE matrices
  u = FEPack.pdes.PDEObject; v = dual(u);
  gradu_gradv = @(muco) (muco * grad2(u)) * grad2(v);
  gradu_vecYv = @(muco) (muco * grad2(u)) * ([0; opts.cutvec(1)] * v);
  vecYu_gradv = @(muco) (muco * ([0; opts.cutvec(1)] * u)) * grad2(v);
  vecYu_vecYv = @(muco) (muco * ([0; opts.cutvec(1)] * u)) * ([0; opts.cutvec(1)] * v);
  u_v = @(rhoco) ((rhoco*u)*v);

  % Edges
  edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
  edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
  edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
  edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

  %% FE matrices
  N_s = meshLineZ.numPoints;
  sVec = meshLineZ.points;
  N_XY = meshXY.numPoints;
  cellXY = meshXY.domain('volumic');
  cutvec = opts.cutvec;
  suffix = opts.suffix;
  omega = opts.omega;
  
  % fid = fopen([nomdossier, 'elapsed_time', suffix], 'w+');
  fid = fopen(['outputs/elapsed_time'], 'a+');

  % parpool
  tic;
  parfor idS = 1:N_s-1
    fprintf('Calcul des matrices éléments finis.\n');
    sVar = sVec(idS, 1);
    mu2Ds  = @(x)  mu3D([x(:, 1), cutvec(1)*x(:, 2), cutvec(2)*x(:, 2) + sVar]); %#ok
    rho2Ds = @(x) rho3D([x(:, 1), cutvec(1)*x(:, 2), cutvec(2)*x(:, 2) + sVar]); %#ok

    FEmat = [];
    FEmat.mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, gradu_gradv(mu2Ds));
    FEmat.mat_gradu_vecYv = FEPack.pdes.Form.intg(cellXY, gradu_vecYv(mu2Ds));
    FEmat.mat_vecYu_gradv = FEPack.pdes.Form.intg(cellXY, vecYu_gradv(mu2Ds));
    FEmat.mat_vecYu_vecYv = FEPack.pdes.Form.intg(cellXY, vecYu_vecYv(mu2Ds));
    FEmat.mat_u_v         = FEPack.pdes.Form.intg(cellXY,        u_v(rho2Ds));

    FEPack.tools.parsave([nomdossier, 'FEmat_', suffix, '_', num2str(idS)], FEmat, true);
  end
  tps = toc;
  fprintf(fid, '%0.5e\t', tps);

  %% Compute Floquet-Bloch transform
  % Homogeneous Dirichlet boundary conditions
  ecs = ((u|edge0x) == 0.0) & ((u|edge1x) == 0.0) &...
        ((u|edge0y) == 0.0) & ((u|edge1y) == 0.0);

  % Surfacic rhs
  B0x = sparse(edge0x.IdPoints(1:end), (1:N0x), 1.0, N_XY, N0x);
  B1x = sparse(edge1x.IdPoints(1:end), (1:N1x), 1.0, N_XY, N1x);
  B0y = sparse(edge0y.IdPoints(2:end-1), (1:N0y-2), 1.0, N_XY, N0y-2);
  B1y = sparse(edge1y.IdPoints(2:end-1), (1:N1y-2), 1.0, N_XY, N1y-2);

  ecs.applyEcs;
  ecs.b = [B0x, B1x, B0y, B1y];
  sidenames = {'0x', '1x', '0y', '1y'};

  % Basis functions on face X = cst and S/Z-line  
  spBX = BCstruct.spBX; NbX = spBX.numBasis; % X = cst
  spBS = BCstruct.spBS; NbS = spBS.numBasis; % S/Z-line
  NbY = NbS * (N0y - 2);

  % Define shear maps
  % /////////////////
  % X = cst
  cut1_0x = opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2);
  cut2_0x = opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2);
  cut1_1x = opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2);
  cut2_1x = opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2);
  %
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([cut1_0x, FEPack.tools.mymod(cut2_0x + s, 0, opts.period), zeros(N0x, 1)]);
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([cut1_1x, FEPack.tools.mymod(cut2_1x + s, 0, opts.period), zeros(N1x, 1)]);
  %
  % S/Z-line
  [IdSX_S, IdSX_X] = ind2sub([NbS, N0y-2], 1:NbY);
  phisY = eye(N0y-2);
  shearmap.fun0y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPack.tools.mymod(s,            0, opts.period) 0, 0], IdSX_S);
  shearmap.fun1y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPack.tools.mymod(s + cutslope, 0, opts.period) 0, 0], IdSX_S);

  % Quadrature for DtN operators: trapezoidal rule
  Nquad = N_s * opts.ratio_quad;
  quadpoints = linspace(0, opts.period, Nquad + 1).';
  quadpoints(end) = [];
  structLoc = meshLineZ.domain('volumic').locateInDomain([quadpoints, zeros(Nquad, 2)]);
  elts = meshLineZ.domain('volumic').elements(structLoc.elements, :);
  coos = structLoc.barycoos;
  quadwgt = (quadpoints(2) - quadpoints(1)) * opts.cutvec(1);
  
  for idFB = 1:opts.numFloquetPoints
    
    % Floquet variable
    FloquetVar = opts.FloquetPoints(idFB);
    fprintf('%d sur %d\n', idFB, opts.numFloquetPoints);
    
    %% Local cell problems
    %  ///////////////////
    fprintf('<strong>Problèmes de cellule locaux.</strong>\n');
    tic;
    parfor idS = 1:N_s-1
      
      % The FE matrix is a linear combination of elementary matrices
      FEmat = load([nomdossier, 'FEmat_', suffix, '_', num2str(idS), '.mat']);

      AA = 1i * FloquetVar * FEmat.mat_vecYu_gradv + FEmat.mat_gradu_gradv...
         - 1i * FloquetVar * FEmat.mat_gradu_vecYv + FloquetVar * FloquetVar * FEmat.mat_vecYu_vecYv...
         - (omega^2) * FEmat.mat_u_v;
      
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
      FEPack.tools.parsave([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(idS)], solcell, true);
    end
    tps = toc;
    fprintf(fid, '%0.5e\t', tps);
    system(['cp ', nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_1.mat ', nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(N_s), '.mat']);
    
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
    parfor idQuad = 1:Nquad

      % Local auxiliary DtN operators
      sVar = quadpoints(idQuad);
      eltsQuad = elts(idQuad, :);
      coosQuad = coos(idQuad, :);
      sc1 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(eltsQuad(1))]);
      sc2 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(eltsQuad(2))]);

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
    
    %% Solve the Riccati equation
    fprintf('<strong>Equation de Riccati.</strong>\n');
    
    flux = @(V) modesFlux(V, E00, E10, F00, F10, orientation, spBX, opts.omega);
    riccatiOpts.tol = 1.0e-2;
    if isfield(opts, 'suffix')
      riccatiOpts.suffix = suffix;
    else
      riccatiOpts.suffix = '';
    end
    riccatiOpts.suffix = [riccatiOpts.suffix, '_', int2str(idFB)];

    % Solve the linearized eigenvalue problem associated to the Riccati equation
    tic;
    [solguide.R, solguide.D] = propagationOperators([E01, E11;  orientation*F01,  orientation*F11], ...
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
    % r_xi_fun = FloquetVar*opts.cutvec(1) + 2*pi*(FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2)).';
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
    
    %% Compute the DtN operator
    U0 = E00 + E10 * solguide.D;
    dU0 = F00 + F10 * solguide.D;

    solguide.Lambda = -BCstruct.BCu' * dU0 + BCstruct.BCdu' * U0;

    %% Save outputs
    save([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)], '-struct', 'solguide', '-v7.3');
  end

  fclose(fid);

end

%% Auxiliary function: The flux function
function flux = modesFlux(V, E00, E10, F00, F10, orientation, spB, omega)
  
  Nb = size(E00, 1);
  Psi  = E00 * V(1:Nb, :) + E10 * V(Nb+1:end, :);
  dPsi = F00 * V(1:Nb, :) + F10 * V(Nb+1:end, :);

  flux = -orientation * imag(diag(Psi' * spB.massmat * dPsi) / omega);

end
