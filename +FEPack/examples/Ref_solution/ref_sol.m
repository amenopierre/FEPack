clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.25i;
problem_setting = 'A'; % 'A' or 'B'
pregenerate_mesh = 0;
struct_mesh = 0;

numCells_pos = 40;
numCells_neg = 40;
BBzmin = -40;
BBzmax =  40;

if strcmpi(problem_setting, 'A')

  % 2D coefficients
  period_posFun = 1;
  period_negFun = 0.5 * sqrt(2);

  mu_pos  = @(x) 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2]);
  rho_pos = @(x) 0.5 + tools.FUN_per_cutoff_cuboid(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu_neg  = @(x) 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)/period_negFun);
  rho_neg = @(x) 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2)/period_negFun);
  
else

  % 2D coefficients
  vecperFun = [-0.5*sqrt(2), 1]; % [-sqrt(2), 1];
  
  mu_pos  = @(x) 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2]);
  rho_pos = @(x) 0.5 + tools.FUN_per_cutoff_cuboid(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu_neg  = @(x) ones(size(x, 1), 1);
  rho_neg = @(x) ones(size(x, 1), 1);

end

% Jump data
alpha_G = 3;
eps_G = 1e-8;
supp_G = -log(eps_G) / alpha_G;
G = @(x) exp(-alpha_G * x(:, 2).^2) .* (abs(x(:, 2)) <= supp_G);
% G = @(x) FEPack.tools.cutoff(x(:, 2), -0.5, 0.5);

%% Mesh
numNodesX = 80;
numNodesZ = 80;

Nz = ceil(numNodesZ * (BBzmax - BBzmin));

if (pregenerate_mesh)

  if (struct_mesh)
    mesh_prefix = 'struct';
  else
    mesh_prefix = 'unstruct';
  end

  % Pick mesh from saved file
  m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodes), '_positive.mat']);
  meshcell = m.mesh;

else

  mesh_pos = meshes.MeshRectangle(struct_mesh, [0  1.0], [BBzmin BBzmax], numNodesX, Nz);
  mesh_neg = meshes.MeshRectangle(struct_mesh, [0 -1.0], [BBzmin BBzmax], numNodesX, Nz);

end


%% Bilinear and linear forms
u = pdes.PDEObject; v = dual(u);
volBilinearIntg = @(muco, rhoco) (muco * grad2(u)) * grad2(v) - (opts.omega^2) * ((rhoco*id(u))*id(v));
volBilinearIntg_pos = volBilinearIntg(mu_pos, rho_pos);
volBilinearIntg_neg = volBilinearIntg(mu_neg, rho_neg);

%% Boundary conditions
basis_functions = 'Lagrange';

% Positive guide
if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domain('xmin'));
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domain('xmax'));
else
  FourierIds = [0 floor(Nz/4)];
  BCstruct_pos.spB0 = spaces.FourierBasis(mesh_pos.domain('xmin'), FourierIds);
  BCstruct_pos.spB1 = spaces.FourierBasis(mesh_pos.domain('xmax'), FourierIds);
end

BCstruct_pos.BCdu = 0.0;
BCstruct_pos.BCu = 1.0;
BCstruct_pos.representation = '';

% Negative guide
if strcmpi(basis_functions, 'Lagrange')
  BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domain('xmin'));
  BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domain('xmax'));
else
  FourierIds = [0 floor(Nz/4)];
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh_neg.domain('xmin'), FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh_neg.domain('xmax'), FourierIds);
end

BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;
BCstruct_neg.representation = '';

%% Compute guide solution
U = PeriodicGuideJumpBVP(1,...
                         volBilinearIntg_pos, mesh_pos, BCstruct_pos, numCells_pos,...
                         volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCells_neg,...
                         G, opts);

%% Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idI = 1:numCells_pos
  trisurf(mesh_pos.triangles, mesh_pos.points(:, 1) + (idI-1),...
                              mesh_pos.points(:, 2), real(U.positive(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  % xlim([-5, 5]);
end

for idI = 1:numCells_neg
  trisurf(mesh_neg.triangles, mesh_neg.points(:, 1) - (idI-1),...
                              mesh_neg.points(:, 2), real(U.negative(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  problem_setting = 'A'; % 'A' or 'B'
  % xlim([-5, 5]);
end
