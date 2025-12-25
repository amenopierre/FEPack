clear; clc;
%%
import FEPack.*

%% Problem-related variables
% //////////////////////////
omega = 20 + 0.25i;
opts.omega = omega;
period = 1;
Lint = 10;
opts.verbose = 0;

% 2D coefficients
vecperfunpos = [cos(3*pi/5), sin(3*pi/5)];
vecperfunneg = [cos(pi/3), sin(pi/3)];

mu2Dpos =  @(x) ones(size(x, 1), 1);
mu2Dneg  = @(x) ones(size(x, 1), 1);

rho2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.4, 0.4]);
rho2Dneg = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperfunneg, [0.5, 0.5], [-0.25, 0.25], [-0.25, 0.25], 0.5, 1);

% mu2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.2, 0.2]);
% rho2Dpos = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
% mu2Dneg  = @(x) ones(size(x, 1), 1); % 1 + 0.5 * cos(2*pi*(x(:, 1) - x(:, 2)*vecperfunneg(1)/vecperfunneg(2))) .* cos(2*pi*x(:, 2)/vecperfunneg(2));
% rho2Dneg = @(x) ones(size(x, 1), 1); % 1 + 0.25 * sin(2*pi*(x(:, 1) - x(:, 2)*vecperfunneg(1)/vecperfunneg(2))) + 0.25 * sin(2*pi*x(:, 2)/vecperfunneg(2));

% Cut matrices
vecperpos = vecperfunpos;
vecperneg = vecperfunneg;

cutvecpos = [1/vecperpos(2); -vecperpos(1)/vecperpos(2)];
cutvecneg = [1/vecperneg(2); -vecperneg(1)/vecperneg(2)];

opts.cutmat_pos = [[1; 0; 0], [0; cutvecpos]];
opts.cutmat_neg = [[1; 0; 0], [0; cutvecneg]];

cutslopePos = cutvecpos(2) / cutvecpos(1);
cutslopeNeg = cutvecneg(2) / cutvecneg(1);

% 3D coefficients
fun3D = @(fun2D, vecper, x) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
mu3Dpos  = @(x) fun3D(mu2Dpos,  vecperpos, x);
mu3Dneg  = @(x) fun3D(mu2Dneg,  vecperneg, x);
rho3Dpos = @(x) fun3D(rho2Dpos, vecperpos, x);
rho3Dneg = @(x) fun3D(rho2Dneg, vecperneg, x);

% Jump data
Gint = @(x) FEPack.tools.cutoff(x(:, 1), -0.5, 0.5);
Gtrs = @(s) ones(size(s, 1), 1);

% (semi-)infinite directions and numbers of cells
semiInfiniteDirection = 1;
infiniteDirection = 2;
numCellsSemiInfinite_pos = 6;
numCellsSemiInfinite_neg = 6;
numCellsInfinite = 5;
numFloquetPoints = 64;

%% Mesh
pregenerate_mesh = 1;
struct_mesh = 1;
numNodes2D = 40;
numNodes3D = 30;

fprintf('Génération maillages\n')
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

% Mesh for the interface
numNodes_int = numNodes2D/4;
mesh2Dint = meshes.MeshSegment('uniform', -Lint, Lint, floor(2*Lint*numNodes_int));
mesh2Dtrs = meshes.MeshSegment('uniform', 0, 1, numNodes2D);

%% Basis functions associated to the interface
% ////////////////////////////////////////////
fprintf('Base associée à l''interface\n')
type_basis_int = 'Lagrange'; % 'sine', 'Lagrange'
type_basis_trs = 'Fourier'; % 'Fourier'

if strcmpi(type_basis_int, 'sine')

  phisInt = @(x, n) sqrt(1 / Lint) * sin(pi * (x(:, 1) + Lint) * n / (2 * Lint)) .* (abs(x(:, 1)) <= Lint);
  numBasisInt = floor(mesh2Dint.numPoints / 2);

elseif strcmpi(type_basis_int, 'Lagrange')

  % error('Il n''y a pas de bug. Je veux juste que tu vérifies le nombre de noeuds (pas trop grand).');
  phisInt = @(x, n) ((x(:, 1) - mesh2Dint.points(n, 1).') ./ (mesh2Dint.points(n+1, 1).' - mesh2Dint.points(n, 1).'))...
                                  .* (x(:, 1) >= mesh2Dint.points(n, 1).' & x(:, 1) < mesh2Dint.points(n+1, 1).')...
                               + ((mesh2Dint.points(n+2, 1).' - x(:, 1)) ./ (mesh2Dint.points(n+2, 1).' - mesh2Dint.points(n+1, 1).'))...
                                  .* (x(:, 1) >= mesh2Dint.points(n+1, 1).' & x(:, 1) < mesh2Dint.points(n+2, 1).');
  numBasisInt = mesh2Dint.numPoints - 2;

else

  error(['Type de fonction de base', type_basis_int, 'non reconnu.']);

end
spBint = spaces.SpectralBasis(mesh2Dint.domain('volumic'), phisInt, numBasisInt);
spBint.computeBasisMatrices;

numBasisTrs = floor(numNodes2D/4);
if strcmpi(type_basis_trs, 'Fourier')

  spBtrs = spaces.FourierBasis(mesh2Dtrs.domain('volumic'), numBasisTrs);

else

  error(['Type de fonction de base', type_basis_trs, 'non reconnu.']);

end

numBasis = spBint.numBasis * spBtrs.numBasis;
[idInt, idTrs] = ind2sub([spBint.numBasis, spBtrs.numBasis], 1:numBasis);

phis3Dpos = @(x, idI) spBtrs.phis(x(:, 3) - x(:, 2) * cutslopePos, idTrs(idI)) .* phisInt(x(:, 2) / cutvecpos(1), idInt(idI));
phis3Dneg = @(x, idI) spBtrs.phis(x(:, 3) - x(:, 2) * cutslopeNeg, idTrs(idI)) .* phisInt(x(:, 2) / cutvecneg(1), idInt(idI));

%%
% phis3Dtest = @(x, idI) spBtrs.phis(x(:, 2) - x(:, 1) * cutslopePos, idTrs(idI)) .* phisInt(x(:, 1) / cutvecpos(1), idInt(idI));
% meshPhi = FEPack.meshes.MeshRectangle(1, [-2 2], [0 1], 40, 40);
% for idI = 1:numBasis
%   trisurf(meshPhi.triangles, meshPhi.points(:, 1), meshPhi.points(:, 2), real(phis3Dtest(meshPhi.points, idI)));
%   view(2);
%   shading interp
%   set(gca, 'DataAspectRatio', [1 1 1])
%   pause;
% end

%% Plot coefficients?
plot_coefficients = false;
if (plot_coefficients)
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  mu2Dglob =  @(x) (x(:, 1) >= 0) .*  mu2Dpos(x) + (x(:, 1) <  0) .* mu2Dneg(x);
  rho2Dglob = @(x) (x(:, 1) >= 0) .* rho2Dpos(x) + (x(:, 1) <  0) .* rho2Dneg(x);

  for idS = 1:(2*numCellsSemiInfinite_pos)
    for idI = 1:(2*numCellsInfinite)
      X = mesh2Dpos.points(:, 1) + (idS - numCellsSemiInfinite_pos - 1);
      Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);
      figure(1);
      trisurf(mesh2Dpos.triangles, X, Y, mu2Dglob([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

      figure(2);
      trisurf(mesh2Dpos.triangles, X, Y, rho2Dglob([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

    end
  end

  figure(1);
  xlim([-numCellsSemiInfinite_pos, numCellsSemiInfinite_pos]);
  ylim([-numCellsInfinite, numCellsInfinite]);

  figure(2);
  xlim([-numCellsSemiInfinite_pos, numCellsSemiInfinite_pos]);
  ylim([-numCellsInfinite, numCellsInfinite]);

  figure(3);
  x = [zeros(256, 1), linspace(-2, 2, 256)', zeros(256, 1)];
  plot(x(:, 2), Gint(x));
  
end

%% Bilinear and linear forms
% //////////////////////////
u = pdes.PDEObject; v = dual(u);
gradu_gradv = @(cutmat, muco) (muco * (cutmat' * grad3(u))) * (cutmat' * grad3(v));
gradu_vec1v = @(cutmat, muco) (muco * (cutmat' * grad3(u))) * (cutmat' * [0; 1; 0] * v);
vec1u_gradv = @(cutmat, muco) (muco * (cutmat' * [0; 1; 0] * u)) * (cutmat' * grad3(v));
vec1u_vec1v = @(cutmat, muco) (muco * (cutmat' * [0; 1; 0] * u)) * (cutmat' * [0; 1; 0] * v);
u_v = @(rhoco) ((rhoco*u)*v);

%% Compute FE elementary matrices
% ///////////////////////////////
fprintf('Matrices EF\n');
% Positive side
tic; mat_gradu_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_gradv(opts.cutmat_pos, mu3Dpos)); toc;
tic; mat_gradu_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_vec1v(opts.cutmat_pos, mu3Dpos)); toc;
tic; mat_vec1u_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_gradv(opts.cutmat_pos, mu3Dpos)); toc;
tic; mat_vec1u_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_vec1v(opts.cutmat_pos, mu3Dpos)); toc;
tic; mat_u_v_pos         = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'),        u_v(rho3Dpos)); toc;

% Negative side
tic; mat_gradu_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_gradv(opts.cutmat_neg, mu3Dneg)); toc;
tic; mat_gradu_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_vec1v(opts.cutmat_neg, mu3Dneg)); toc;
tic; mat_vec1u_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_gradv(opts.cutmat_neg, mu3Dneg)); toc;
tic; mat_vec1u_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_vec1v(opts.cutmat_neg, mu3Dneg)); toc;
tic; mat_u_v_neg         = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'),        u_v(rho3Dneg)); toc;

%% Boundary conditions
% ////////////////////
basis_functions = 'Lagrange';
if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmin'));
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmax'));
else
  FourierIds = [0, numNodes3D/4, numNodes3D/4];
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
  FourierIds = [0, numNodes3D/4, numNodes3D/4];
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh3Dneg.domain('xmin'), FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh3Dneg.domain('xmax'), FourierIds);
end
BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;
BCstruct_neg.representation = 'projection';


%% Floquet-Bloch transform of the positive half-space solution
% ////////////////////////////////////////////////////////////
FloquetPoints_pos = linspace(-pi/period, pi/period, numFloquetPoints).';
opts.computeSol = false;
sol_pos_data = cell(numFloquetPoints, 1);
TFBlambdaPos = cell(numFloquetPoints, 1);
parfor idFB = 1:numFloquetPoints
  fprintf('%d sur %d\n', idFB, numFloquetPoints);
  
  FloquetVar = FloquetPoints_pos(idFB);

  % The FE matrix is a linear combination of the elementary pieces
  AApos = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vec1u_gradv_pos - 1i * FloquetVar * mat_gradu_vec1v_pos + FloquetVar * FloquetVar * mat_vec1u_vec1v_pos - (omega^2) * mat_u_v_pos;

  % Compute the Floquet-Bloch transform of the solution
  [~, sol_pos_data{idFB}.E0, sol_pos_data{idFB}.E1, sol_pos_data{idFB}.R, sol_pos_data{idFB}.D, ~, TFBlambdaPos{idFB}] = PeriodicHalfGuideBVP(mesh3Dpos, +1, semiInfiniteDirection, AApos, BCstruct_pos, numCellsSemiInfinite_pos, opts);
end

%% Floquet-Bloch transform of the negative half-space solution
% ////////////////////////////////////////////////////////////
FloquetPoints_neg = linspace(-pi/period, pi/period, numFloquetPoints).';
opts.computeSol = false;
sol_neg_data = cell(numFloquetPoints, 1);
TFBlambdaNeg = cell(numFloquetPoints, 1);
parfor idFB = 1:numFloquetPoints
  fprintf('%d sur %d\n', idFB, numFloquetPoints);
  
  FloquetVar = FloquetPoints_neg(idFB);

  % The FE matrix is a linear combination of the elementary pieces
  AAneg = mat_gradu_gradv_neg + 1i * FloquetVar * mat_vec1u_gradv_neg - 1i * FloquetVar * mat_gradu_vec1v_neg + FloquetVar * FloquetVar * mat_vec1u_vec1v_neg - (omega^2) * mat_u_v_neg;

  % Compute the Floquet-Bloch transform of the solution
  [~, sol_neg_data{idFB}.E0, sol_neg_data{idFB}.E1, sol_neg_data{idFB}.R, sol_neg_data{idFB}.D, ~, TFBlambdaNeg{idFB}] = PeriodicHalfGuideBVP(mesh3Dneg, +1, semiInfiniteDirection, AAneg, BCstruct_neg, numCellsSemiInfinite_neg, opts);
end

%% Positive half-space DtN
% ////////////////////////
TFBphiPosVec = zeros(mesh3Dpos.domain('xmin').numPoints, numBasis, numFloquetPoints);
for idI = 1:numBasis
  TFBphiPosVec(:, idI, :) = FEPack.tools.BlochTransform(mesh3Dpos.points(mesh3Dpos.domain('xmin').IdPoints, :),...
                                           FloquetPoints_pos, @(x) phis3Dpos(x, idI),...
                                           2, period, 1000);
end

% %%
% meshDom = mesh3Dpos.restrictToDomain(mesh3Dpos.domain('xmin'));
% for idI = 1:numBasis
%   trisurf(meshDom.triangles, meshDom.points(:, 2), meshDom.points(:, 3), real(TFBphiPosVec(:, idI, 2)));
%   view(2);
%   shading interp
%   set(gca, 'DataAspectRatio', [1 1 1])
%   pause;
% end
% %%

TFBphiPos = cell(numFloquetPoints, 1);
lambda_pos = zeros(numBasis);
wpos = (2*pi / period) / (numFloquetPoints - 1);

for idFB = 1:numFloquetPoints
  TFBphiPos{idFB} = BCstruct_pos.spB0.FE_to_spectral * TFBphiPosVec(:, :, idFB);
  TFBlambdaPos{idFB} = BCstruct_pos.spB0.massmat * TFBlambdaPos{idFB};

  lambda_pos = lambda_pos + wpos * TFBphiPos{idFB}' * TFBlambdaPos{idFB} * TFBphiPos{idFB};
end

lambda_pos = lambda_pos / cutvecpos(1);

%% Negative half-space DtN
% ////////////////////////
TFBphiNegVec = zeros(mesh3Dneg.domain('xmin').numPoints, numBasis, numFloquetPoints);
for idI = 1:numBasis
  TFBphiNegVec(:, idI, :) = FEPack.tools.BlochTransform(mesh3Dneg.points(mesh3Dneg.domain('xmin').IdPoints, :),...
                                           FloquetPoints_neg, @(x) phis3Dneg(x, idI),...
                                           2, period, 1000);
end

TFBphiNeg = cell(numFloquetPoints, 1);
lambda_neg = zeros(numBasis);
wneg = (2*pi / period) / (numFloquetPoints - 1);

for idFB = 1:numFloquetPoints
  TFBphiNeg{idFB} = BCstruct_neg.spB0.FE_to_spectral * TFBphiNegVec(:, :, idFB);
  TFBlambdaNeg{idFB} = BCstruct_neg.spB0.massmat * TFBlambdaNeg{idFB};

  lambda_neg = lambda_neg + wneg * TFBphiNeg{idFB}' * TFBlambdaNeg{idFB} * TFBphiNeg{idFB};
end

lambda_neg = lambda_neg / cutvecneg(1);

%% Solve the integral equation on the interface
% /////////////////////////////////////////////
numNodesEF = mesh2Dint.numPoints * mesh2Dtrs.numPoints;
[idIntEF, idTrsEF] = ind2sub([mesh2Dint.numPoints, mesh2Dtrs.numPoints], 1:numNodesEF);

projmat = spBint.projmat(idInt, idIntEF) .* spBtrs.projmat(idTrs, idTrsEF);
Gvec = Gint(mesh2Dint.points) * Gtrs(mesh2Dtrs.points).';
Gvec = Gvec(:);
mat_G_v = projmat * Gvec;

% -------------------------------------------------------------------- %
% The minus sign comes from the definition of the Lambda when they are %
% the DtN operators (see PeriodicHalfGuideBVP.m)                       %
% -------------------------------------------------------------------- %
trace_solution = -(lambda_pos + lambda_neg) \ mat_G_v;             %
% -------------------------------------------------------------------- %

%% Deduce the FB transform of the solution
% ////////////////////////////////////////
% FB transform of the solution's trace with period_pos
TFB_solution_pos = cell(numFloquetPoints, 1);
for idFB = 1:numFloquetPoints

  TFB_solution_pos{idFB} = zeros(mesh3Dpos.numPoints, numCellsSemiInfinite_pos);
  R0Phi = TFBphiPos{idFB} * trace_solution;
  R1Phi = sol_pos_data{idFB}.D * R0Phi;

  for idCell = 0:numCellsSemiInfinite_pos-1
    % Compute the solution in the current cell
    TFB_solution_pos{idFB}(:, idCell + 1) = sol_pos_data{idFB}.E0 * R0Phi + sol_pos_data{idFB}.E1 * R1Phi;

    % Update
    R0Phi = sol_pos_data{idFB}.R * R0Phi;
    R1Phi = sol_pos_data{idFB}.D * R0Phi;
  end

end

% FB transform of the solution's trace with period_neg
TFB_solution_neg = cell(numFloquetPoints, 1);
for idFB = 1:numFloquetPoints

  TFB_solution_neg{idFB} = zeros(mesh3Dneg.numPoints, numCellsSemiInfinite_neg);
  R0Phi = TFBphiNeg{idFB} * trace_solution;
  R1Phi = sol_neg_data{idFB}.D * R0Phi;

  for idCell = 0:numCellsSemiInfinite_neg-1
    % Compute the solution in the current cell
    TFB_solution_neg{idFB}(:, idCell + 1) = sol_neg_data{idFB}.E0 * R0Phi + sol_neg_data{idFB}.E1 * R1Phi;

    % Update
    R0Phi = sol_neg_data{idFB}.R * R0Phi;
    R1Phi = sol_neg_data{idFB}.D * R0Phi;
  end

end

%% Take the inverse Floquet transform: positive side
% //////////////////////////////////////////////////
numCells = [1 1 1];
numCells(infiniteDirection) = 2*numCellsInfinite;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_pos;
Nu = prod(numCells);
[I1, I2, I3] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2; I3]; % 3-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period) ./ (numFloquetPoints - 1)); % 1-by-1
U3D.positive = zeros(mesh3Dpos.numPoints, Nu);

for idFB = 1:numFloquetPoints
  FloquetVar = FloquetPoints_pos(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dpos.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFB_solution_pos{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.positive = U3D.positive + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.positive = U3D.positive * sqrt(period / (2*pi));

%% Take the inverse Floquet transform: negative side
% //////////////////////////////////////////////////
numCells = [1 1 1];
numCells(infiniteDirection) = 2*numCellsInfinite;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_neg;
Nu = prod(numCells);
[I1, I2, I3] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2; I3]; % 3-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period) ./ (numFloquetPoints - 1)); % 1-by-1
U3D.negative = zeros(mesh3Dneg.numPoints, Nu);

for idFB = 1:numFloquetPoints
  FloquetVar = FloquetPoints_neg(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dneg.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFB_solution_neg{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.negative = U3D.negative + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.negative = U3D.negative * sqrt(period / (2*pi));

%% Take the trace: Positive side
% //////////////////////////////
N2Dpos = mesh2Dpos.numPoints;
U2D.positive = zeros(N2Dpos, size(U3D.positive, 2));
dom = mesh3Dpos.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_pos*(idI-1)+1):(numCellsSemiInfinite_pos*idI);
  X = mesh2Dpos.points(:, 1);
  Y = mesh2Dpos.points(:, 2); % ones(mesh2Dpos.numPoints, 1);
  Z = FEPack.tools.mymod(cutslopePos * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.positive(:, IcellY) = reshape(sum(reshape(coos .* U3D.positive(elts, IcellY), dom.dimension+1, []), 1), N2Dpos, []);
end

%% Take the trace: Negative side
% //////////////////////////////
N2Dneg = mesh2Dneg.numPoints;
U2D.negative = zeros(N2Dneg, size(U3D.negative, 2));
dom = mesh3Dneg.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_neg*(idI-1)+1):(numCellsSemiInfinite_neg*idI);
  X = mesh2Dneg.points(:, 1);
  Y = mesh2Dneg.points(:, 2); % ones(mesh2Dneg.numPoints, 1);
  Z = FEPack.tools.mymod(cutslopeNeg * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.negative(:, IcellY) = reshape(sum(reshape(coos .* U3D.negative(elts, IcellY), dom.dimension+1, []), 1), N2Dneg, []);
end

% profile viewer
% profile OFF

%% Plot U
% ///////
mafig = figure;%(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% if (compareU)
%   subplot(1, 2, 1);
% end
for idS = 1:numCellsSemiInfinite_pos
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
    X = mesh2Dpos.points(:, 1) + (idS - 1);
    Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);

    trisurf(mesh2Dpos.triangles, X, Y, real(U2D.positive(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
    % caxis([-0.02, 0.02]);
  end
end
%%
for idS = 1:numCellsSemiInfinite_neg
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
    X = mesh2Dneg.points(:, 1) - (idS - 1);
    Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
    trisurf(mesh2Dneg.triangles, X, Y, real(U2D.negative(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
  end
end

xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
ylim([-numCellsInfinite, numCellsInfinite]);

savefig(mafig, 'outputs/solution', 'compact');
