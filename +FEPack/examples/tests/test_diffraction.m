clear; clc;
%%
import FEPack.*

mu_neg  = @(x) ones(size(x, 1), 1);
rho_neg = @(x) ones(size(x, 1), 1);
mu_pos  = @(x) ones(size(x, 1), 1);
rho_pos = @(x) 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], [0; 1], [0.5, 0.5], [-0.2, 0.2]);

opts.omega = 4 + 0.1i;
opts.computeSol = true;
opts.solBasis = true;

numCells_pos = 10;
numCells_neg = 10;

% Incident field
angle_diff = pi/3;
vec_diff = (opts.omega^2) * rho_neg(0) * [cos(angle_diff); sin(angle_diff)];
uinc  = @(x) exp(1i * x(:, 1:2) * vec_diff);
dx_uinc = @(x) 1i * vec_diff(1) * uinc(x);

% Mesh
numNodesX = 4;
numNodesY = numNodesX;

mesh_pos = meshes.MeshRectangle(1, [0  1], [0, 1], numNodesX, numNodesY);
cellXY = mesh_pos.domain('volumic');
edge0x = mesh_pos.domain('xmin');
edge1x = mesh_pos.domain('xmax');

u = pdes.PDEObject; 
v = dual(u);

mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, (mu_pos  * (grad2(u))) * (grad2(v)));
mat_gradu_vectv = FEPack.pdes.Form.intg(cellXY, (mu_pos  * (grad2(u))) * ([0; 1] * v));
mat_vectu_gradv = FEPack.pdes.Form.intg(cellXY, (mu_pos  * ([0; 1] * u)) * (grad2(v)));
mat_vectu_vectv = FEPack.pdes.Form.intg(cellXY, (mu_pos  * ([0; 1] * u)) * ([0; 1] * v));
mat_u_v         = FEPack.pdes.Form.intg(cellXY, (rho_pos * u) * v);

AApos = mat_gradu_gradv...
      + 1i * vec_diff(2)     * mat_vectu_gradv...
      - 1i * vec_diff(2)     * mat_gradu_vectv...
      + (vec_diff(2) * vec_diff(2)) * mat_vectu_vectv...
      + mat_u_v;

BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(edge0x);
BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(edge1x);
BCstruct_pos.BCdu = 0.0;
BCstruct_pos.BCu = 1.0;
BCstruct_pos.representation = '';

[~, E0pos, E1pos, Rpos, Dpos, BCstruct_pos, Lambda_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 1, AApos, BCstruct_pos, numCells_pos, opts);

% Negative side
mesh_neg = meshes.MeshRectangle(0, [0 -1], [0, 1], numNodesX, numNodesY);
cellXY = mesh_neg.domain('volumic');
edge0x = mesh_neg.domain('xmin');
edge1x = mesh_neg.domain('xmax');

u = pdes.PDEObject; 
v = dual(u);

mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, (mu_neg  * (grad2(u))) * (grad2(v)));
mat_gradu_vectv = FEPack.pdes.Form.intg(cellXY, (mu_neg  * (grad2(u))) * ([0; 1] * v));
mat_vectu_gradv = FEPack.pdes.Form.intg(cellXY, (mu_neg  * ([0; 1] * u)) * (grad2(v)));
mat_vectu_vectv = FEPack.pdes.Form.intg(cellXY, (mu_neg  * ([0; 1] * u)) * ([0; 1] * v));
mat_u_v         = FEPack.pdes.Form.intg(cellXY, (rho_neg * u) * v);

AAneg = mat_gradu_gradv...
      + 1i * vec_diff(2)     * mat_vectu_gradv...
      - 1i * vec_diff(2)     * mat_gradu_vectv...
      + (vec_diff(2) * vec_diff(2)) * mat_vectu_vectv...
      + mat_u_v;

BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(edge0x);
BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(edge1x);
BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;
BCstruct_neg.representation = '';

[~, E0neg, E1neg, Rneg, Dneg, BCstruct_neg, Lambda_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 1, AAneg, BCstruct_neg, numCells_neg, opts);

% Interface equation
xOnInt = [zeros(edge0x.numPoints, 1), mesh_pos.points(edge0x.IdPoints, 2)];
Lambda_pos = -Lambda_pos;
Lambda_neg = -Lambda_neg;

GG =             - BCstruct_neg.spB0.FE_to_spectral * dx_uinc(xOnInt) +...
      Lambda_neg * BCstruct_neg.spB0.FE_to_spectral *    uinc(xOnInt);

psiTrans = (Lambda_pos + Lambda_neg) \ GG;
psiDiff  = psiTrans - BCstruct_pos.spB0.FE_to_spectral * uinc(xOnInt);


%% Construct the solution in the whole domain
% //////////////////////////////////////////
% psi = BCstruct_pos.spB0.FE_to_spectral * ones(edge0x.numPoints, 1);

% Positive side
U.positive = zeros(mesh_pos.numPoints, numCells_pos);
R0Phi = psiTrans;
R1Phi = Dpos * R0Phi;

for idCell = 0:numCells_pos-1
  % Compute the solution in the current cell
  U.positive(:, idCell + 1) = E0pos * R0Phi + E1pos * R1Phi;

  % Update the coefficients
  R0Phi = Rpos * R0Phi;
  R1Phi = Dpos * R0Phi;
end

% Negative side
U.negative = zeros(mesh_neg.numPoints, numCells_neg);
R0Phi = psiDiff;
R1Phi = Dneg * R0Phi;

for idCell = 0:numCells_neg-1
  X = [mesh_neg.points(:, 1) - idCell, mesh_neg.points(:, 2)];

  % Compute the solution in the current cell
  % U.negative(:, idCell + 1) = E0neg * R0Phi + E1neg * R1Phi + uinc(X);
  U.negative(:, idCell + 1) = E0neg * R0Phi + E1neg * R1Phi;

  % Update the coefficients
  R0Phi = Rneg * R0Phi;
  R1Phi = Dneg * R0Phi;
end

% Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idI = 1:numCells_pos
  trisurf(mesh_pos.triangles, mesh_pos.points(:, 1) + (idI-1),...
                              mesh_pos.points(:, 2), real(U.positive(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
end

for idI = 1:numCells_neg
  trisurf(mesh_neg.triangles, mesh_neg.points(:, 1) - (idI-1),...
                              mesh_neg.points(:, 2), real(U.negative(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
end

% xlim([-5, 5]);
