% clear; clc;
import FEPack.*

% profile ON
opts.omega = 20 + 0.25i;
problem_setting = 'B'; % 'A' or 'B'
pregenerate_mesh = 1;
struct_mesh = 1;

numCellsXpos = 10;
numCellsXneg = 10;
sizeCellZ = 1.0;
Zorigin = -sizeCellZ * 0.5 * numCellsZ;

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
% alpha_G = 100;
% eps_G = 1e-8;
% supp_G = -log(eps_G) / alpha_G;
% G = @(x) exp(-alpha_G * x(:, 2).^2) .* (abs(x(:, 2)) <= supp_G);
G = @(x) 100 * FEPack.tools.cutoff(x(:, 2), -0.5, 0.5);

%% Initialization
fprintf('A. Initialisation\n');

% Mesh
numNodesX = 500;
numNodesZ = 500;

Nz = ceil(numNodesZ * sizeCellZ);

if (pregenerate_mesh)

  if (struct_mesh)
    mesh_prefix = 'struct';
  else
    mesh_prefix = 'unstruct';
  end

  % Pick mesh from saved file
  m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodesX), '_positive.mat']);
  mesh_pos = m.mesh;

  m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodesX), '_negative_X.mat']);
  mesh_neg = m.mesh;

else

  mesh_pos = meshes.MeshRectangle(struct_mesh, [0  1.0], [0 sizeCellZ], numNodesX, Nz);
  mesh_neg = meshes.MeshRectangle(struct_mesh, [0 -1.0], [0 sizeCellZ], numNodesX, Nz);

end

basis_function.name = 'Fourier';
basis_function.FourierIds = (-numNodesX:numNodesX);

% Save the workspace
save([cheminDonnees, '/inputs.mat'], '-v7.3');