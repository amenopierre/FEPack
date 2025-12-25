% clear; clc;
%%

% problem_setting = entrees{1};
% period          = entrees{2};
% function transmission_quasi2D_par(problem_setting, period, numNodes, test_omega, test_case, test_per, test_G, compute_3D_sol, compare_w_refsol, delete_sol, suffix)
  import FEPack.*
  profile OFF

  nargin = 0;
  
  %% Default variables
  if (nargin < 1)
    problem_setting = 'A'; % 'A' or 'B'
  end
  if (nargin < 2)
    period = 1;
  end
  if (nargin < 3)
    numNodes = 10;
  end
  if (nargin < 4)
    test_omega = 'omega1';
  end
  if (nargin < 5)
    test_case = 'constant'; % 'constant' or 'rational' or 'irrational' or 'homogeneisation'
  end
  if (nargin < 6)
    test_per = 1;
  end
  if (nargin < 7)
    test_G = 'G1';
  end
  if (nargin < 8)
    compute_3D_sol = false;
  end
  if (nargin < 9)
    compare_w_refsol = false;
  end
  if (nargin < 10)
    delete_sol = false;
  end
  if (nargin < 11)
    suffix = '';
  end

  % opts.nomdossier = ['last_simus/solution_', problem_setting, '_', test_case, '_', test_omega, '_per', num2str(test_per), '_', test_G, suffix, '/']; % /!\ Pas oublier le / à la fin
  opts.nomdossier = ['/UMA/tmp/amenoagbadji/solution_', problem_setting, '_', test_case, '_', test_omega, '_per', num2str(test_per), '_', test_G, suffix, '/']; % /!\ Pas oublier le / à la fin
  nomdossier = opts.nomdossier;

  u = FEPack.pdes.PDEObject; v = dual(u);
  disp(dy(u));
  error();

  system(['mkdir ', nomdossier]);
  
  if ~(compare_w_refsol)
    warning(['Fichiers dans ', nomdossier, ' supprimés.']);
    system(['rm ', nomdossier, '*']);
  end

  %% Problem-related variables
  opts.period = period;
  opts.verbose = 0;
  
  if strcmpi(test_omega, 'omega1')
    opts.omega = 8 + 0.25i;
  elseif strcmpi(test_omega, 'omega2')
    opts.omega = 20 + 0.25i;
  elseif strcmpi(test_omega, 'omega3')
    opts.omega = 20 + 0.05i;
  elseif strcmpi(test_omega, 'omega4')
    opts.omega = 1 + 0.25i;
  end

  mu2Dpos =  @(x) ones(size(x, 1), 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1);

  if strcmpi(problem_setting, 'A')

    % 2D coefficients
    period_posFun = 1;

    if strcmpi(test_case, 'rational')
      period_negFun = 1;
    else
      period_negFun = sqrt(2);
    end

    % echelle = 0.1;
    if strcmpi(test_case, 'constant')
      rho2Dpos = @(x) ones(size(x, 1), 1);
      rho2Dneg = @(x) 1 * ones(size(x, 1), 1);
    else
      % funperiod = opts.period;
      % opts.period = 1;

      rho2DposCell = @(x) 0.5 + perCutoffCircle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.4, 0.4]);
      rho2DnegCell = @(x) 0.5 + perCutoffCuboid(x, [1; 0], [0; period_negFun], [0.5, 0.5], [-0.25, 0.25], [-0.25, 0.25], 0.5, 1);
      
      % rho_eff_pos = integral2(@(X, Y) reshape(rho2DposCell([X(:), Y(:)]), size(X)), 0, 1, 0, period_posFun) / period_posFun;
      % rho_eff_neg = integral2(@(X, Y) reshape(rho2DnegCell([X(:), Y(:)]), size(X)), 0, 1, 0, period_negFun) / period_negFun;

      rho2Dpos = @(x) rho2DposCell(x / opts.period);
      rho2Dneg = @(x) rho2DnegCell(x / opts.period);
    end

    % Cut vector
    period_pos = period_posFun;
    period_neg = test_per * period_negFun;

    opts.cutvec = [1/period_pos; 1/period_neg];

    % 3D coefficients
    mu3Dpos  = @(x) mu2Dpos( [x(:, 1), period_pos * x(:, 2)]);
    mu3Dneg  = @(x) mu2Dneg( [x(:, 1), period_neg * x(:, 3)]);
    rho3Dpos = @(x) rho2Dpos([x(:, 1), period_pos * x(:, 2)]);
    rho3Dneg = @(x) rho2Dneg([x(:, 1), period_neg * x(:, 3)]);
    
  else

    % 2D coefficients
    if strcmpi(test_case, 'rational')
      vecperFun = [1, 1];
    else
      vecperFun = [cos(3*pi/5), sin(3*pi/5)];
    end
    
    if strcmpi(test_case, 'constant')
      rho2Dpos = @(x) ones(size(x, 1), 1);
      rho2Dneg = @(x) 2 * ones(size(x, 1), 1);
    else
      rho2Dpos = @(x) 0.5 + perCutoffCircle(x/opts.period, [1; 0], vecperFun, [0.5, 0.5], [-0.4, 0.4]);
      rho2Dneg = @(x) ones(size(x, 1), 1);
    end

    % Cut vector
    vecper = test_per * [-vecperFun(1), vecperFun(2)];
    opts.cutvec = [1/vecper(2); -vecper(1)/vecper(2)];
    
    % 3D coefficients
    fun3D = @(fun2D, x) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
    mu3Dpos = @(x) fun3D(mu2Dpos, x);
    mu3Dneg = @(x) fun3D(mu2Dneg, x);
    rho3Dpos = @(x) fun3D(rho2Dpos, x);
    rho3Dneg = @(x) fun3D(rho2Dneg, x);

  end

  % Cut matrix and cut slope
  opts.cutmat = [[1; 0; 0], [0; opts.cutvec]];

  % Jump data
  G = @(x) ones(size(x, 1), 1); % FEPack.tools.cutoff(x, -0.5, 0.5);

  if strcmpi(test_G, 'G1')
    G3D = @(x) G(x(:, 1)/opts.cutvec(1));
  else
    G3D = @(x) G(x(:, 1)/opts.cutvec(1)) .* exp(2i*pi*(x(:, 2) - x(:, 1)*opts.cutvec(2)/opts.cutvec(1)));
  end

  % Numbers of cells
  numCellsSemiInfinite_pos = ceil(1/opts.period);
  numCellsSemiInfinite_neg = ceil(1/opts.period);
  opts.numCellsInfinite = ceil(1*opts.cutvec(1)/opts.period);
  opts.numFloquetPoints = 1;% ceil(50/opts.period);
  opts.numFloquetPoints = max(opts.numFloquetPoints, ceil(opts.numFloquetPoints/opts.cutvec(1)));
  opts.FloquetPoints = 0;%linspace(-pi/period, pi/period, opts.numFloquetPoints);

  %% Meshes
  struct_mesh = 1;

  numNodesX = numNodes; % numNodesX = max(10, ceil(numNodesX*opts.period))
  numNodesY = numNodes; % numNodesY = max(10, ceil(numNodesY/opts.cutvec(1)))
  numNodesZ = numNodes;
  opts.ratio_quad = 1;

  meshXYpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period/opts.cutvec(1)], numNodesX, numNodesY);
  meshXYneg = FEPack.meshes.MeshRectangle(struct_mesh, [0 -opts.period], [0  opts.period/opts.cutvec(1)], numNodesX, numNodesY);
  meshZXpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesZ, numNodesX);
  meshZXneg = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0 -opts.period], numNodesZ, numNodesX);
  meshYZ    = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesY, numNodesZ);
  meshLineZ = FEPack.meshes.MeshSegment('uniform', 0, opts.period, numNodesZ);

  % %% Plot the coefficients
  % mafig = figure;
  % set(groot,'defaultAxesTickLabelInterpreter','latex');
  % set(groot,'defaulttextinterpreter','latex');
  % set(groot,'defaultLegendInterpreter','latex');

  % for idY = 1:2*opts.numCellsInfinite
  %   for idX = 1:numCellsSemiInfinite_pos
  %     X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
  %     Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

  %     trisurf(meshXYpos.triangles, X, Y, rho2Dpos([X, Y]));
  %     hold on;
  %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %   end
  % end

  % for idY = 1:2*opts.numCellsInfinite
  %   for idX = 1:numCellsSemiInfinite_neg
  %     X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
  %     Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

  %     trisurf(meshXYneg.triangles, X, Y, rho2Dneg([X, Y]));
  %     hold on;
  %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %   end
  % end

  % savefig(mafig, [nomdossier, 'fig_coeff_2D'], 'compact');
  % close(mafig);

  %% Basis functions and boundary conditions
  % Positive side
  basis_functions_X = 'Fourier';
  basis_functions_Y = 'Lagrange';
  if strcmpi(basis_functions_X, 'Lagrange')
    BCstruct.pos.spBX = FEPack.spaces.PeriodicLagrangeBasis(meshYZ.domain('volumic')); % X = cst
  else
    nBY = max(1, floor(numNodesY/2));
    nBZ = max(1, floor(numNodesZ/2));
    BCstruct.pos.spBX = FEPack.spaces.FourierBasis(meshYZ.domain('volumic'), [nBY, nBZ]);
  end

  if strcmpi(basis_functions_Y, 'Lagrange')
    % Particular Lagrange for basis functons on Y = cst
    ecs_fun = @(u, dom) ((((u|dom.mesh.domain('xmin')) - (u|dom.mesh.domain('xmax'))) == 0.0) &...      % Periodic on Z
                          ((u|dom.mesh.domain('ymin')) == 0.0) &...
                          ((u|dom.mesh.domain('ymax')) == 0.0)); % Dirichlet on X
    
    BCstruct.pos.spBY = FEPack.spaces.EssentialBCLagrangeBasis(meshZXpos.domain('volumic'), ecs_fun); % Y = cst
  else
    nBZ = max(1, floor(numNodesZ/4));
    nBX = max(1, floor(numNodesX/4));
    
    % SBidsZ = (-nBZ:nBZ).';
    % SBidsX = (1:nBX).';
    % numBZX = (2*nBZ + 1) * nBX;
    
    % BCstruct.pos.spBY = FEPack.spaces.SpectralBasis(meshZXpos.domain('volumic'), @(P, n) spDirichletPeriodic(P, n, SBidsZ, SBidsX, [opts.period, opts.period]), numBZX);
    % BCstruct.pos.spBY.computeBasisMatrices;
    BCstruct.pos.spBY = FEPack.spaces.FourierBasis(meshZXpos.domain('volumic'), [nBZ, nBX]); % Y = cst
  end

  % Other parameters
  BCstruct.pos.BCdu = 0.0;
  BCstruct.pos.BCu = 1.0;
  % BCstruct.pos.representation = 'projection';

  % Negative side
  BCstruct.neg = BCstruct.pos;

  if strcmpi(basis_functions_Y, 'Lagrange')
    % Particular Lagrange for basis functons on Y = cst
    ecs_fun = @(u, dom) ((((u|dom.mesh.domain('xmin')) - (u|dom.mesh.domain('xmax'))) == 0.0) &...      % Periodic on Z
                          ((u|dom.mesh.domain('ymin')) == 0.0) & ((u|dom.mesh.domain('ymax')) == 0.0)); % Dirichlet on X

    BCstruct.neg.spBY = FEPack.spaces.EssentialBCLagrangeBasis(meshZXneg.domain('volumic'), ecs_fun); % Y = cst
  else
    nBZ = max(1, floor(numNodesZ/4));
    nBX = max(1, floor(numNodesX/4));

    % SBidsZ = (-nBZ:nBZ).';
    % SBidsX = (1:nBX).';
    % numBZX = (2*nBZ + 1) * nBX;
    
    % BCstruct.neg.spBY = FEPack.spaces.SpectralBasis(meshZXneg.domain('volumic'), @(P, n) spDirichletPeriodic(P, n, SBidsZ, SBidsX, [opts.period, opts.period]), numBZX);
    % BCstruct.neg.spBY.computeBasisMatrices;
    BCstruct.neg.spBY = FEPack.spaces.FourierBasis(meshZXneg.domain('volumic'), [nBZ, nBX]); % Y = cst
  end

  %% Compute DtN operators
  opts.suffix = 'pos';
  quasi2DHalfGuideDirichlet_par(+1, meshXYpos, meshLineZ, mu3Dpos, rho3Dpos, BCstruct.pos, opts);
  system(['rm ', nomdossier, 'FEmat_*']);
  %%
  opts.suffix = 'neg';
  quasi2DHalfGuideDirichlet_par(-1, meshXYneg, meshLineZ, mu3Dneg, rho3Dneg, BCstruct.neg, opts);
  system(['rm ', nomdossier, 'FEmat_*']);

  %% Interface equation
  fprintf('<strong>Equation d''interface.</strong>\n');
  FloquetPoints = opts.FloquetPoints;
  numFloquetPoints = opts.numFloquetPoints;
  pointsYZ = meshYZ.points;
  spBX_FE_to_spectral = BCstruct.pos.spBX.FE_to_spectral;
  period = opts.period;

  % Nv = (size(spBX_FE_to_spectral, 1) + 1)/2;
  FourierIdsX = BCstruct.pos.spBX.FourierIds.X(:).'; dx = length(FourierIdsX);
  FourierIdsY = BCstruct.pos.spBX.FourierIds.Y(:).'; dy = length(FourierIdsY);
  FourierIdsZ = BCstruct.pos.spBX.FourierIds.Z(:).'; dz = length(FourierIdsZ);
  Iconst = sub2ind([dx dy dz], find(FourierIdsX==0), find(FourierIdsY==0), find(FourierIdsZ==0));

  parfor idFB = 1:numFloquetPoints
    fprintf('%d sur %d\n', idFB, numFloquetPoints);
    FloquetVar = FloquetPoints(idFB);

    % DtN operators
    solguidePos = load([nomdossier, 'half_guide_sol_pos_Floquet_', num2str(idFB)]);
    solguideNeg = load([nomdossier, 'half_guide_sol_neg_Floquet_', num2str(idFB)]);
    
    Lambda_pos = solguidePos.Lambda;
    Lambda_neg = solguideNeg.Lambda;

    % Rhs
    % jumpData_FB = @(x) FEPack.tools.BlochTransform(x, FloquetVar, G3D, 1, period, 1000);
    % GG = spBX_FE_to_spectral * jumpData_FB(pointsYZ);
    
    GG = zeros(size(spBX_FE_to_spectral, 1), 1); GG(Iconst) = 1;% spBX_FE_to_spectral * ones(meshYZ.numPoints, 1);
      
    % The minus sign comes from the definition of the Lambda
    soltrace = [];
    % soltrace.vec = -(Lambda_pos + Lambda_neg) \ GG;
    soltrace.vec = GG;

    % Save the trace of the solution
    FEPack.tools.parsave([nomdossier, 'sol_trace_Floquet_', num2str(idFB)], soltrace, true);
  end

  %% Construct solution
  opts.compute3Dsolution = compute_3D_sol;
  % numNodesX = 10; numNodesY = 10; numNodesZ = 10;
  opts.meshYZ = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesY, numNodesZ);
  mesh3DXYpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesX, numNodesY);
  mesh3DZXpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesZ, numNodesX);
  mesh3DXYneg = FEPack.meshes.MeshRectangle(struct_mesh, [0 -opts.period], [0  opts.period], numNodesX, numNodesY);
  mesh3DZXneg = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0 -opts.period], numNodesZ, numNodesX);

  opts.suffix = 'pos';
  opts.meshXY = mesh3DXYpos;
  opts.meshZX = mesh3DZXpos;
  opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
  construct_solution_par(+1, meshXYpos, meshLineZ, BCstruct.pos, opts);

  opts.suffix = 'neg';
  opts.meshXY = mesh3DXYneg;
  opts.meshZX = mesh3DZXneg;
  opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
  construct_solution_par(-1, meshXYneg, meshLineZ, BCstruct.neg, opts);

  %% Plot solution
  if strcmpi(test_case, 'homogeneisation')
    close all;
    numCellsInfinite = opts.numCellsInfinite;

    %
    % Compute limit solution
    % //////////////////////
    xibound = 30;
    xi = linspace(-xibound, xibound, floor(100 * (2 * xibound))).';
    Delta_xi = xi(2) - xi(1);

    % Fourier transform of interface data
    ypts = linspace(-0.5, 0.5, ceil(10*xibound)).';
    Delta_y = ypts(2) - ypts(1);
    FT_G = Delta_y * exp(-1i * xi * ypts.') * G(ypts);
    
    % //////////////////// %
    % Configuration A only %
    % //////////////////// %
    rho_eff_pos = integral2(@(X, Y) reshape(rho2DposCell([X(:), Y(:)]), size(X)), 0, 1, 0, period_posFun) / period_posFun;
    rho_eff_neg = integral2(@(X, Y) reshape(rho2DnegCell([X(:), Y(:)]), size(X)), 0, 1, 0, period_negFun) / period_negFun;
    %
    coe_pos = sqrt(xi.^2 - rho_eff_pos * opts.omega^2);
    coe_neg = sqrt(xi.^2 - rho_eff_neg * opts.omega^2);
    %
    solrefPos = @(X, Y) (Delta_xi/(2*pi)) * exp(-X * coe_pos.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));
    solrefNeg = @(X, Y) (Delta_xi/(2*pi)) * exp(+X * coe_neg.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));
    %
    %
    %
    mafig = figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_pos
        X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
        Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

        trisurf(meshXYpos.triangles, X, Y, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end
    %
    %
    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_neg
        X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
        Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

        trisurf(meshXYneg.triangles, X, Y, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end

    xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
    ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
    savefig(mafig, [nomdossier, 'fig_solution_2D'], 'compact');
    close(mafig);


    mafig = figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_pos
        X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
        Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

        solErreur = sol.Usol(:, idX) - solrefPos(X, Y);
        trisurf(meshXYpos.triangles, X, Y, real(solErreur));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end
    %
    %
    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_neg
        X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
        Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

        solErreur = sol.Usol(:, idX) - solrefNeg(X, Y);
        trisurf(meshXYneg.triangles, X, Y, real(solErreur));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end

    xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
    ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
    savefig(mafig, [nomdossier, 'fig_erreur_2D'], 'compact');
    close(mafig);

    % meshPosPoints = [];
    % meshPosTriangles = [];
    % tri = meshXYpos.triangles;
    % UsolPos = [];
    % UsolPosRef = [];

    % parfor idY = 1:2*numCellsInfinite
    %   sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);
      
    %   for idX = 1:numCellsSemiInfinite_pos
    %     fprintf('%d/%d et %d/%d\n', idX, numCellsSemiInfinite_pos, idY, 2*numCellsInfinite);
    %     idpt = (idY-1) * numCellsSemiInfinite_pos + (idX-1);

    %     X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
    %     Y = meshXYpos.points(:, 2) + (idY - numCellsInfinite - 1) * opts.period / opts.cutvec(1);

    %     meshPosPoints = [meshPosPoints; [X, Y]];
    %     meshPosTriangles = [meshPosTriangles; tri+idpt*meshXYpos.numPoints];

    %     UsolPos = [UsolPos; real(sol.Usol(:, idX))]
    %     UsolPosRef = [UsolPosRef; real(solrefPos(X, Y))];
    %   end
    % end

    % meshNegPoints = [];
    % meshNegTriangles = [];
    % tri = meshXYneg.triangles;
    % UsolNeg = [];
    % UsolNegRef = [];

    % parfor idY = 1:2*numCellsInfinite
    %   sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);
      
    %   for idX = 1:numCellsSemiInfinite_neg
    %     fprintf('%d/%d et %d/%d\n', idX, numCellsSemiInfinite_neg, idY, 2*numCellsInfinite);
    %     idpt = (idY-1) * numCellsSemiInfinite_neg + (idX-1);

    %     X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
    %     Y = meshXYneg.points(:, 2) + (idY - numCellsInfinite - 1) * opts.period / opts.cutvec(1);

    %     meshNegPoints = [meshNegPoints; [X, Y]];
    %     meshNegTriangles = [meshNegTriangles; tri+idpt*meshXYneg.numPoints];

    %     UsolNeg = [UsolNeg; real(sol.Usol(:, idX))]
    %     UsolNegRef = [UsolNegRef; real(solrefNeg(X, Y))];
    %   end
    % end

    % mafig = figure;
    % set(groot,'defaultAxesTickLabelInterpreter','latex');
    % set(groot,'defaulttextinterpreter','latex');
    % set(groot,'defaultLegendInterpreter','latex');

    % trisurf(meshPosTriangles, meshPosPoints(:, 1), meshPosPoints(:, 2), UsolPos); hold on;
    % trisurf(meshNegTriangles, meshNegPoints(:, 1), meshNegPoints(:, 2), UsolNeg);

    % savefig(mafig, [nomdossier, 'fig_solution_2D'], 'compact');
    % close(mafig);
    % %
    % mafig = figure;
    % set(groot,'defaultAxesTickLabelInterpreter','latex');
    % set(groot,'defaulttextinterpreter','latex');
    % set(groot,'defaultLegendInterpreter','latex');

    % trisurf(meshPosTriangles, meshPosPoints(:, 1), meshPosPoints(:, 2), UsolPos-UsolPosRef); hold on;
    % trisurf(meshNegTriangles, meshNegPoints(:, 1), meshNegPoints(:, 2), UsolNeg-UsolNegRef);

    % savefig(mafig, [nomdossier, 'fig_erreur_2D'], 'compact');
    % close(mafig);

  elseif (strcmpi(test_case, 'constant') & compare_w_refsol)
    % Compute reference solution
    % //////////////////////
    xibound = 30;
    xi = 0; % linspace(-xibound, xibound, floor(500 * (2 * xibound))).';
    % Delta_xi = xi(2) - xi(1);

    % Fourier transform of interface data
    % ypts = linspace(-0.5, 0.5, ceil(10*xibound)).';
    % Delta_y = ypts(2) - ypts(1);
    FT_G = 1;% Delta_y * exp(-1i * xi * ypts.') * G(ypts);
    
    % //////////////////// %
    % Configuration A only %
    % //////////////////// %
    rho_const_pos = rho2Dpos([0, 0]);
    rho_const_neg = rho2Dneg([0, 0]);
    %
    xi = 0;
    coe_pos = sqrt(xi.^2 - rho_const_pos * opts.omega^2);
    coe_neg = sqrt(xi.^2 - rho_const_neg * opts.omega^2);
    %
    % solrefPos = @(X, Y) (Delta_xi/(2*pi)) * exp(-X * coe_pos.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));
    % solrefNeg = @(X, Y) (Delta_xi/(2*pi)) * exp(+X * coe_neg.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));
    solrefPos = @(X, Y) exp(-X * coe_pos); % (Delta_xi/(2*pi)) * exp(-X * coe_pos.' + 1i * Y * xi.') * FT_G;
    solrefNeg = @(X, Y) exp(+X * coe_neg); % (Delta_xi/(2*pi)) * exp(+X * coe_neg.' + 1i * Y * xi.') * FT_G;

    u = FEPack.pdes.PDEObject; v = dual(u);
    mat_gradu_gradv_pos = FEPack.pdes.Form.intg(meshXYpos.domain('volumic'), (dy(u) * dy(v)));
    mat_u_v_pos         = FEPack.pdes.Form.intg(meshXYpos.domain('volumic'), u * v);
    mat_gradu_gradv_neg = FEPack.pdes.Form.intg(meshXYneg.domain('volumic'), grad2(u) * grad2(v));
    mat_u_v_neg         = FEPack.pdes.Form.intg(meshXYneg.domain('volumic'), u * v);
    
    ecartL2pos = 0; ecartH1pos = 0;
    normUL2pos = 0; normUH1pos = 0;
    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);
      for idX = 1:numCellsSemiInfinite_pos
        X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
        Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);
        % X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
        % Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);
        ecart = sol.Usol(:, idX) - solrefPos(X, Y);

        ecartL2pos = ecartL2pos + (ecart' * mat_u_v_pos * ecart);
        normUL2pos = normUL2pos + (solrefPos(X, Y)' * mat_u_v_pos * solrefPos(X, Y));

        ecartH1pos = ecartH1pos + (ecart' * mat_gradu_gradv_pos * ecart);
        % ecartH1pos = ecartH1pos + (ecart' * mat_u_v_pos * ecart);
        normUH1pos = normUH1pos + (solrefPos(X, Y)' * mat_gradu_gradv_pos * solrefPos(X, Y));
      end
    end
    fid = fopen([nomdossier, 'erreurs_pos.txt'], 'a+');
    fprintf(fid, '%d\t%0.5d\t%0.5d\t%0.5d\t%0.5d\n', numNodesX, ecartL2pos, normUL2pos, ecartH1pos, normUH1pos);
    fclose(fid);

    ecartL2neg = 0; ecartH1neg = 0;
    normUL2neg = 0; normUH1neg = 0;
    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);
      for idX = 1:numCellsSemiInfinite_neg
        X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
        Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);
        ecart = sol.Usol(:, idX) - solrefNeg(X, Y);

        ecartL2neg = ecartL2neg + (ecart' * mat_u_v_neg * ecart);
        normUL2neg = normUL2neg + (solrefNeg(X, Y)' * mat_u_v_neg * solrefNeg(X, Y));

        ecartH1neg = ecartH1neg + (ecart' * mat_gradu_gradv_neg * ecart);
        normUH1neg = normUH1neg + (solrefNeg(X, Y)' * mat_gradu_gradv_neg * solrefNeg(X, Y));
      end
    end
    fid = fopen([nomdossier, 'erreurs_neg.txt'], 'a+');
    fprintf(fid, '%d\t%0.5d\t%0.5d\t%0.5d\t%0.5d\n', numNodesX, ecartL2neg, normUL2neg, ecartH1neg, normUH1neg);
    fclose(fid);

  else
    mafig = figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_pos
        X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
        Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

        trisurf(meshXYpos.triangles, X, Y, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end
    %
    %
    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_neg
        X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
        Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

        trisurf(meshXYneg.triangles, X, Y, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end

    xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
    ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
    % ylim([1/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
    % caxis([-0.25, 0.25]);
    savefig(mafig, [nomdossier, 'fig_solution_2D'], 'compact');
  end

  system(['rm ', nomdossier, 'half_guide_*']);
  system(['rm ', nomdossier, 'local_cell_sol_*']);
  system(['rm ', nomdossier, 'sol_trace_Floquet_*']);
  if (delete_sol)
    system(['rm ', nomdossier, 'solution_*']);
  end


  %% Plot 3D solution
  if isfield(opts, 'compute3Dsolution') && (opts.compute3Dsolution)
    
    % X constant
    for idX = 1:numCellsSemiInfinite_pos      
      mafig3D_X = figure;
      set(groot,'defaultAxesTickLabelInterpreter','latex');
      set(groot,'defaulttextinterpreter','latex');
      set(groot,'defaultLegendInterpreter','latex');
      
      for idY = 1:2*opts.numCellsInfinite
        sol = load([nomdossier, 'solution3D_Xcst_pos_Y_', num2str(idY)]);
        
        Y = opts.meshYZ.points(:, 1) + (idY - opts.numCellsInfinite - 1) * opts.period;
        Z = opts.meshYZ.points(:, 2);

        trisurf(opts.meshYZ.triangles, Y, Z, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
      end
      
      savefig(mafig3D_X, [nomdossier, 'fig_solution_3D_pos_X_', int2str(idX)]);
      close(mafig3D_X);
    end

    for idX = 1:numCellsSemiInfinite_neg      
      mafig3D_X = figure;
      set(groot,'defaultAxesTickLabelInterpreter','latex');
      set(groot,'defaulttextinterpreter','latex');
      set(groot,'defaultLegendInterpreter','latex');
      
      for idY = 1:2*opts.numCellsInfinite
        sol = load([nomdossier, 'solution3D_Xcst_neg_Y_', num2str(idY)]);
        
        Y = opts.meshYZ.points(:, 1) + (idY - opts.numCellsInfinite - 1) * opts.period;
        Z = opts.meshYZ.points(:, 2);

        trisurf(opts.meshYZ.triangles, Y, Z, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
      end
      
      savefig(mafig3D_X, [nomdossier, 'fig_solution_3D_neg_X_', int2str(idX)]);
      close(mafig3D_X);
    end

    % Y constant
    for idY = 1:2*opts.numCellsInfinite
      mafig3D_Y = figure;
      set(groot,'defaultAxesTickLabelInterpreter','latex');
      set(groot,'defaulttextinterpreter','latex');
      set(groot,'defaultLegendInterpreter','latex');

      % Côté positif  
      sol = load([nomdossier, 'solution3D_Ycst_pos_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_pos
        Z = mesh3DZXpos.points(:, 1);
        X = mesh3DZXpos.points(:, 2) + (idX - 1) * opts.period;

        trisurf(mesh3DZXpos.triangles, Z, X, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end

      % Côté négatif  
      sol = load([nomdossier, 'solution3D_Ycst_neg_Y_', num2str(idY)]);

      for idX = 1:numCellsSemiInfinite_pos
        Z = mesh3DZXneg.points(:, 1);
        X = mesh3DZXneg.points(:, 2) - (idX - 1) * opts.period;

        trisurf(mesh3DZXneg.triangles, Z, X, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
      
      savefig(mafig3D_Y, [nomdossier, 'fig_solution_3D_Y_', int2str(idY)]);
      close(mafig3D_Y);
    end
    

    % Z constant
    mafig3D_Z = figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    
    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution3D_Zcst_pos_Y_', num2str(idY)]);
    
      for idX = 1:numCellsSemiInfinite_pos
        X = mesh3DXYpos.points(:, 1) + (idX - 1) * opts.period;
        Y = mesh3DXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period;
    
        trisurf(mesh3DXYpos.triangles, X, Y, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end
    %
    %
    for idY = 1:2*opts.numCellsInfinite
      sol = load([nomdossier, 'solution3D_Zcst_neg_Y_', num2str(idY)]);
    
      for idX = 1:numCellsSemiInfinite_neg
        X = mesh3DXYneg.points(:, 1) - (idX - 1) * opts.period;
        Y = mesh3DXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period;
    
        trisurf(mesh3DXYneg.triangles, X, Y, real(sol.Usol(:, idX)));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end

    savefig(mafig3D_Z, [nomdossier, 'fig_solution_3D_Z']);
    close(mafig3D_Z);
    %
    system(['rm ', nomdossier, 'solution3D_*']);
  end

% end


function val = spDirichletPeriodic(P, n, IdsX, IdsY, periods)
    IdsX = IdsX(:).'; dx = length(IdsX);
    IdsY = IdsY(:).'; dy = length(IdsY);

    [Ix, Iy] = ind2sub([dx, dy], n(:)');

    val = exp(2i*pi*P(:, 1)*IdsX(Ix)/periods(1)) .* sin(pi*P(:, 2)*IdsY(Iy)/periods(2));
end