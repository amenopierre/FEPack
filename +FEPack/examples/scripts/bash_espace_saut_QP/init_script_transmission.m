function init_script_transmission(numFloquetPoints, numNodes, cheminDonnees)
  
  if (nargin < 1)
    numFloquetPoints = 50;
  end

  if (nargin < 2)
    numNodes = 10;
  end

  %%
  import FEPack.*
  % profile ON

  %% Problem-related variables
  opts.omega = 8 + 0.25i;
  period = 1;
  opts.verbose = 0;
  problem_setting = 'B'; % 'A' or 'B'

  if strcmpi(problem_setting, 'A')

    % 2D coefficients
    period_posFun = 1;
    period_negFun = sqrt(2);

    mu2Dpos  = @(x) ones(size(x, 1), 1); % 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2]);
    rho2Dpos = @(x) ones(size(x, 1), 1); % 0.5 + tools.FUN_per_cutoff_cuboid(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
    mu2Dneg  = @(x) ones(size(x, 1), 1); % 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)/period_negFun);
    rho2Dneg = @(x) ones(size(x, 1), 1); % 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2)/period_negFun);

    % Cut vector
    period_pos = period_posFun;
    period_neg = period_negFun;

    cutvec = [period_pos; period_neg];

    % 3D coefficients
    mu3Dpos = @(x) mu2Dpos([x(:, 1), period_pos * x(:, 2)]);
    rho3Dpos = @(x) rho2Dpos([x(:, 1), period_pos * x(:, 2)]);
    mu3Dneg = @(x) mu2Dneg([x(:, 1), period_neg * x(:, 3)]);
    rho3Dneg = @(x) rho2Dneg([x(:, 1), period_neg * x(:, 3)]);
    
  else

    % 2D coefficients
    vecperFun = [-sqrt(2), 1]; % [-sqrt(2), 1];
    
    mu2Dpos = @(x) 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2]);
    rho2Dpos = @(x) 0.5 + tools.FUN_per_cutoff_cuboid(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
    mu2Dneg  = @(x) ones(size(x, 1), 1);
    rho2Dneg = @(x) ones(size(x, 1), 1);

    % Cut vector
    vecper = vecperFun;
    cutvec = [1/vecper(2); -vecper(1)/vecper(2)];
    
    % 3D coefficients
    fun3D = @(fun2D, x) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
    mu3Dpos = @(x) fun3D(mu2Dpos, x);
    mu3Dneg = @(x) fun3D(mu2Dneg, x);
    rho3Dpos = @(x) fun3D(rho2Dpos, x);
    rho3Dneg = @(x) fun3D(rho2Dneg, x);

  end

  % Cut matrix and cut slope
  opts.cutmat = [[1; 0; 0], [0; cutvec]];
  cutslope = cutvec(2) / cutvec(1);

  % Jump data
  % alpha_G = 3;
  % eps_G = 1e-8;
  % supp_G = -log(eps_G) / alpha_G;
  % G = @(x) exp(-alpha_G * x(:, 2).^2) .* (abs(x(:, 2)) <= supp_G);
  G = @(x) 10 * FEPack.tools.cutoff(x(:, 2), -0.5, 0.5);
  G3D = @(x) G([zeros(size(x, 1), 1), x(:, 2)/cutvec(1), zeros(size(x, 1), 1)]) .* exp(2i*pi*(x(:, 3) - x(:, 2)*cutvec(2)/cutvec(1)));

  % (semi-)infinite directions and numbers of cells
  semiInfiniteDirection = 1;
  infiniteDirection = 2;
  numCellsSemiInfinite_pos = 10;
  numCellsSemiInfinite_neg = 10;
  numCellsInfinite = 10;

  %% Mesh
  pregenerate_mesh = 1;
  struct_mesh = 1;
  numNodes2D = 50;
  numNodes3D = numNodes;

  if pregenerate_mesh

    if (struct_mesh)
      mesh_prefix = 'struct';
    else
      mesh_prefix = 'unstruct';
    end

    % Pick mesh from saved file
    m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodes2D), '_positive.mat']);
    mesh2Dpos = m.mesh;

    m = load(['pregenMeshes/3D/', mesh_prefix, '_mesh_3D_', int2str(numNodes3D), '_positive.mat']);
    mesh3Dpos = m.mesh;

    m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodes2D), '_negative_X.mat']);
    mesh2Dneg = m.mesh;

    m = load(['pregenMeshes/3D/', mesh_prefix, '_mesh_3D_', int2str(numNodes3D), '_negative_X.mat']);
    mesh3Dneg = m.mesh;

  else

    mesh2Dpos = meshes.MeshRectangle(struct_mesh, [0 1], [0 1], numNodes2D, numNodes2D);
    mesh3Dpos = meshes.MeshCuboid(struct_mesh, [0 1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

    mesh2Dneg = meshes.MeshRectangle(struct_mesh, [0 -1], [0 1], numNodes2D, numNodes2D);
    mesh3Dneg = meshes.MeshCuboid(struct_mesh, [0 -1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

  end

  %% Bilinear and linear forms
  u = pdes.PDEObject; v = dual(u);
  gradu_gradv = @(muco) (muco * (opts.cutmat' * grad3(u))) * (opts.cutmat' * grad3(v));
  gradu_vec1v = @(muco) (muco * (opts.cutmat' * grad3(u))) * (opts.cutmat' * [0; 1; 0] * v);
  vec1u_gradv = @(muco) (muco * (opts.cutmat' * [0; 1; 0] * u)) * (opts.cutmat' * grad3(v));
  vec1u_vec1v = @(muco) (muco * (opts.cutmat' * [0; 1; 0] * u)) * (opts.cutmat' * [0; 1; 0] * v);
  u_v = @(rhoco) ((rhoco*u)*v);

  %% Compute FE elementary matrices
  % Positive side
  tic; mat_gradu_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_gradv(mu3Dpos)); toc;
  tic; mat_gradu_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_vec1v(mu3Dpos)); toc;
  tic; mat_vec1u_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_gradv(mu3Dpos)); toc;
  tic; mat_vec1u_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_vec1v(mu3Dpos)); toc;
  tic; mat_u_v_pos         = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'),        u_v(rho3Dpos)); toc;

  % Negative side
  tic; mat_gradu_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_gradv(mu3Dneg)); toc;
  tic; mat_gradu_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_vec1v(mu3Dneg)); toc;
  tic; mat_vec1u_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_gradv(mu3Dneg)); toc;
  tic; mat_vec1u_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_vec1v(mu3Dneg)); toc;
  tic; mat_u_v_neg         = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'),        u_v(rho3Dneg)); toc;

  %% Boundary conditions
  basis_functions = 'Lagrange';
  if strcmpi(basis_functions, 'Lagrange')
    BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmin'));
    BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmax'));
  else
    FourierIds = [0, floor(numNodes3D/2), 0];
    BCstruct_pos.spB0 = spaces.FourierBasis(mesh3Dpos.domain('xmin'), FourierIds);
    BCstruct_pos.spB1 = spaces.FourierBasis(mesh3Dpos.domain('xmax'), FourierIds);
  end
  BCstruct_pos.BCdu = 0.0;
  BCstruct_pos.BCu = 1.0;
  BCstruct_pos.representation = 'projection';

  if strcmpi(basis_functions, 'Lagrange')
    BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmin'));
    BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmax'));
  else
    FourierIds = [0, floor(numNodes3D/2), 0];
    BCstruct_neg.spB0 = spaces.FourierBasis(mesh3Dneg.domain('xmin'), FourierIds);
    BCstruct_neg.spB1 = spaces.FourierBasis(mesh3Dneg.domain('xmax'), FourierIds);
  end
  BCstruct_neg.BCdu = 0.0;
  BCstruct_neg.BCu = 1.0;
  BCstruct_neg.representation = 'projection';

  % Floquet points
  FloquetPoints = linspace(-pi/period, pi/period, numFloquetPoints);

  %% Save the workspace
  save([cheminDonnees, '/inputs_', int2str(numNodes), '.mat'], '-v7.3');
  
end