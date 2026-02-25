clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.25i;
problem_setting = 'A'; % 'A' or 'B'
pregenerate_mesh = 0;
struct_mesh = 0;

numCellsXpos = 6;
numCellsXneg = 6;
numCellsZ = 5;
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
alpha_G = 3;
eps_G = 1e-8;
supp_G = -log(eps_G) / alpha_G;
G = @(x) exp(-alpha_G * x(:, 2).^2) .* (abs(x(:, 2)) <= supp_G);
% G = @(x) FEPack.tools.cutoff(x(:, 2), -0.5, 0.5);

%% Initialization
fprintf('A. Initialisation\n');

% Mesh
numNodesX = 20;
numNodesZ = 20;

Nz = ceil(numNodesZ * sizeCellZ);

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

  mesh_pos = meshes.MeshRectangle(struct_mesh, [0  1.0], [0 sizeCellZ], numNodesX, Nz);
  mesh_neg = meshes.MeshRectangle(struct_mesh, [0 -1.0], [0 sizeCellZ], numNodesX, Nz);

end

%% Trace of guide solution
fprintf('B. Problèmes de demi-guide\n');
referenceHalfGuide(mesh_pos, +1, opts.omega, mu_pos, rho_pos, numCellsZ, sizeCellZ, Zorigin, 'pos');
referenceHalfGuide(mesh_neg, -1, opts.omega, mu_neg, rho_neg, numCellsZ, sizeCellZ, Zorigin, 'neg');

fprintf('C. Transmission à l''interface\n');
solGuidePos = load('outputs/half_guide_solution_pos');
solGuideNeg = load('outputs/half_guide_solution_neg');

% Right-hand side
numBasisX = size(solGuidePos.Lambda, 1);
Sigma0x = mesh_pos.domain('xmin'); N0x = Sigma0x.numPoints;
u = FEPack.pdes.PDEObject; v = dual(u);
GG = sparse(numBasisX, 1);

MM = FEPack.pdes.Form.intg(mesh_pos.domain('volumic'), u * v);
MM = MM(Sigma0x.IdPoints, Sigma0x.IdPoints);

for idZ = 1:numCellsZ
  II = (1+(idZ-1)*(N0x-1):1+idZ*(N0x-1))';
  cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;

  VV = G([mesh_pos.points(Sigma0x.IdPoints, 1), mesh_pos.points(Sigma0x.IdPoints, 2) + cellZorigin]);
  GG = GG + sparse(II, 1, VV, numBasisX, 1);
end

% Solve linear system
solphi = (solGuidePos.Lambda + solGuideNeg.Lambda) \ GG;

%% Compute guide solution
Usol.positive = cell(numCellsXpos, numCellsZ);
Usol.negative = cell(numCellsXneg, numCellsZ);

fprintf('D. Construction de la solution\n');
for idZ = 1:numCellsZ

  fprintf('%3d sur %3d\n', idZ, numCellsZ);
  
  % Solution in the positive guide
  R0phi = solphi;
  R1phi = solGuidePos.S * solphi;
  SC = load(['outputs/local_cell_sol_pos_', int2str(idZ)]);
  idBasis = 1+(idZ-1)*(N0x-1):1+idZ*(N0x-1);

  for idX = 1:numCellsXpos
    Usol.positive{idX, idZ} = (SC.E0x * solGuidePos.B0(idBasis, :) +...
                               SC.E0z * solGuidePos.traceE0{idZ} +...
                               SC.E1z * solGuidePos.traceE0{idZ+1}) * R0phi...
                               ...
                            + (SC.E1x * solGuidePos.B1(idBasis, :) +...
                               SC.E0z * solGuidePos.traceE1{idZ} +...
                               SC.E1z * solGuidePos.traceE1{idZ+1}) * R1phi;
    
    % Update the coefficients
    R0phi = solGuidePos.P * R0phi;
    R1phi = solGuidePos.S * R0phi;
  end

  % Solution in the negative guide
  R0phi = solphi;
  R1phi = solGuideNeg.S * solphi;
  SC = load(['outputs/local_cell_sol_neg_', int2str(idZ)]);
  idBasis = 1+(idZ-1)*(N0x-1):1+idZ*(N0x-1);

  for idX = 1:numCellsXneg
    Usol.negative{idX, idZ} = (SC.E0x * solGuideNeg.B0(idBasis, :) +...
                               SC.E0z * solGuideNeg.traceE0{idZ} +...
                               SC.E1z * solGuideNeg.traceE0{idZ+1}) * R0phi...
                               ...
                            + (SC.E1x * solGuideNeg.B1(idBasis, :) +...
                               SC.E0z * solGuideNeg.traceE1{idZ} +...
                               SC.E1z * solGuideNeg.traceE1{idZ+1}) * R1phi;
    
    % Update the coefficients
    R0phi = solGuideNeg.P * R0phi;
    R1phi = solGuideNeg.S * R0phi;
  end

end

save('outputs/guide_solution');

%% Plot Guide solution
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idX = 1:numCellsXpos
  for idZ = 1:numCellsZ
    fprintf('%d/%d et %d/%d\n', idX, numCellsXpos, idZ, numCellsZ);
    cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;

    trisurf(mesh_pos.triangles, mesh_pos.points(:, 1) + (idX - 1),...
                                mesh_pos.points(:, 2) + cellZorigin, real(Usol.positive{idX, idZ}));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

for idX = 1:numCellsXneg
  for idZ = 1:numCellsZ
    fprintf('%d/%d et %d/%d\n', idX, numCellsXneg, idZ, numCellsZ);
    cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;

    trisurf(mesh_neg.triangles, mesh_neg.points(:, 1) - (idX - 1),...
                                mesh_neg.points(:, 2) + cellZorigin, real(Usol.negative{idX, idZ}));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end
