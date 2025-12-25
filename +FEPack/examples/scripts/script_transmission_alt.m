clear; clc;
%%
import FEPack.*
% profile ON

%% Problem-related variables
opts.omega = 8 + 0.5i;
opts.verbose = 0;
opts.computeSol = true;
opts.solBasis = true;

period_pos = 1;
period_neg = 0.5*sqrt(2);
Lint = 10;

mu2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], [0; period_pos], [0.5, 0.5], [-0.2, 0.2]);
rho2Dpos = @(x) 0.5 + perCutoffCuboid(x, [1; 0], [0; period_pos], [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
mu2Dneg  = @(x) 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)/period_neg);
rho2Dneg = @(x) 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2)/period_neg);

% Jump data
alpha_G = 3;
eps_G = 1e-8;
supp_G = -log(eps_G) / alpha_G;
Gint = @(x) exp(-alpha_G * x(:, 1).^2) .* (abs(x(:, 1)) <= supp_G);
% Gint = @(x) FEPack.tools.cutoff(x(:, 1), -0.1, 0.1);

% (semi-)infinite directions and numbers of cells
semiInfiniteDirection = 1;
infiniteDirection = 2;
numCellsSemiInfinite_pos = 6;
numCellsSemiInfinite_neg = 6;
numCellsInfinite_pos = 6;
numCellsInfinite_neg = 9;
numFloquetPoints_pos = 50;
numFloquetPoints_neg = 50;

numFloquetPoints_pos = floor(numFloquetPoints_pos / period_pos);
numFloquetPoints_neg = floor(numFloquetPoints_neg / period_neg);


%% Mesh
struct_mesh = 1;
numNodes_int = 20;
numNodes_pos = 20;
numNodes_neg = 20;

mesh2Dpos = meshes.MeshRectangle(struct_mesh, [0 1], [0 period_pos], numNodes_pos, floor(period_pos * numNodes_pos));
mesh2Dneg = meshes.MeshRectangle(struct_mesh, [0 -1], [0 period_neg], numNodes_neg, floor(period_neg * numNodes_neg));
mesh2Dint = meshes.MeshSegment('uniform', -Lint, Lint, floor(2*Lint*numNodes_int));

%% Basis functions associated to the interface
type_basis = 'Lagrange'; % 'sine', 'Lagrange'

if strcmpi(type_basis, 'sine')

  int_basis_functions = @(x, n) sqrt(1 / Lint) * sin(pi * (x(:, 1) + Lint) * n / (2 * Lint)) .* (abs(x(:, 1)) <= Lint);
  numBasisInt = 20;

elseif strcmpi(type_basis, 'Lagrange')

  % error('Il n''y a pas de bug. Je veux juste que tu vérifies le nombre de noeuds (pas trop grand).');
  int_basis_functions = @(x, n) ((x(:, 1) - mesh2Dint.points(n, 1).') ./ (mesh2Dint.points(n+1, 1).' - mesh2Dint.points(n, 1).'))...
                                  .* (x(:, 1) >= mesh2Dint.points(n, 1).' & x(:, 1) < mesh2Dint.points(n+1, 1).')...
                               + ((mesh2Dint.points(n+2, 1).' - x(:, 1)) ./ (mesh2Dint.points(n+2, 1).' - mesh2Dint.points(n+1, 1).'))...
                                  .* (x(:, 1) >= mesh2Dint.points(n+1, 1).' & x(:, 1) < mesh2Dint.points(n+2, 1).');
  numBasisInt = mesh2Dint.numPoints - 2;

else

  error(['Type de fonction de base', type_basis, 'non reconnu.']);

end

spBint = spaces.SpectralBasis(mesh2Dint.domain('volumic'), int_basis_functions, numBasisInt);
spBint.computeBasisMatrices;

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
  plot(x(:, 2), G(x));
end

%% Bilinear and linear forms
u = pdes.PDEObject; v = dual(u);
gradu_gradv = @(muco) (muco * grad2(u)) * grad2(v);
gradu_vecZv = @(muco) (muco * grad2(u)) * ([0; 1] * v);
vecZu_gradv = @(muco) (muco * ([0; 1] * u)) * grad2(v);
vecZu_vecZv = @(muco) (muco * ([0; 1] * u)) * ([0; 1] * v);
u_v = @(rhoco) (rhoco * u) * v;

%% Compute FE elementary matrices
% Positive side
mat_gradu_gradv_pos = FEPack.pdes.Form.intg(mesh2Dpos.domain('volumic'), gradu_gradv(mu2Dpos));
mat_gradu_vecZv_pos = FEPack.pdes.Form.intg(mesh2Dpos.domain('volumic'), gradu_vecZv(mu2Dpos));
mat_vecZu_gradv_pos = FEPack.pdes.Form.intg(mesh2Dpos.domain('volumic'), vecZu_gradv(mu2Dpos));
mat_vecZu_vecZv_pos = FEPack.pdes.Form.intg(mesh2Dpos.domain('volumic'), vecZu_vecZv(mu2Dpos));
mat_u_v_pos         = FEPack.pdes.Form.intg(mesh2Dpos.domain('volumic'),        u_v(rho2Dpos));

% Negative side
mat_gradu_gradv_neg = FEPack.pdes.Form.intg(mesh2Dneg.domain('volumic'), gradu_gradv(mu2Dneg));
mat_gradu_vecZv_neg = FEPack.pdes.Form.intg(mesh2Dneg.domain('volumic'), gradu_vecZv(mu2Dneg));
mat_vecZu_gradv_neg = FEPack.pdes.Form.intg(mesh2Dneg.domain('volumic'), vecZu_gradv(mu2Dneg));
mat_vecZu_vecZv_neg = FEPack.pdes.Form.intg(mesh2Dneg.domain('volumic'), vecZu_vecZv(mu2Dneg));
mat_u_v_neg         = FEPack.pdes.Form.intg(mesh2Dneg.domain('volumic'),        u_v(rho2Dneg));

%% Boundary conditions
basis_functions = 'Lagrange';
if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmin'));
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmax'));
else
  FourierIds = [0, numNodes2D/2, 0];
  BCstruct_pos.spB0 = spaces.FourierBasis(mesh2Dpos.domain('xmin'), FourierIds);
  BCstruct_pos.spB1 = spaces.FourierBasis(mesh2Dpos.domain('xmax'), FourierIds);
end
BCstruct_pos.BCdu = 0.0;
BCstruct_pos.BCu = 1.0;
BCstruct_pos.representation = 'weak evaluation';

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmin'));
  BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmax'));
else
  FourierIds = [0, numNodes2D/2, 0];
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh2Dneg.domain('xmin'), FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh2Dneg.domain('xmax'), FourierIds);
end
BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;
BCstruct_neg.representation = 'weak evaluation';

%% Floquet-Bloch transform of the positive half-space solution
FloquetPoints_pos = linspace(-pi/period_pos, pi/period_pos, numFloquetPoints_pos)';
sol_pos_data = cell(numFloquetPoints_pos, 1);
TFBlambdaPos = cell(numFloquetPoints_pos, 1);

for idFB = 1:numFloquetPoints_pos
  fprintf('%d sur %d\n', idFB, numFloquetPoints_pos);
  
  FloquetVar = FloquetPoints_pos(idFB);

  % The FE matrix is a linear combination of the elementary pieces
  AApos = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vecZu_gradv_pos - 1i * FloquetVar * mat_gradu_vecZv_pos + FloquetVar * FloquetVar * mat_vecZu_vecZv_pos  - (opts.omega^2) * mat_u_v_pos;

  % Compute the Floquet-Bloch transform of the solution
  [~, sol_pos_data{idFB}.E0, sol_pos_data{idFB}.E1, sol_pos_data{idFB}.R, sol_pos_data{idFB}.D, BCstruct_pos, TFBlambdaPos{idFB}] = PeriodicHalfGuideBVP(mesh2Dpos, 1, semiInfiniteDirection, AApos, BCstruct_pos, numCellsSemiInfinite_pos, opts);
end

%% Floquet-Bloch transform of the negative half-space solution
FloquetPoints_neg = linspace(-pi/period_neg, pi/period_neg, numFloquetPoints_neg)';
sol_neg_data = cell(numFloquetPoints_neg, 1);
TFBlambdaNeg = cell(numFloquetPoints_neg, 1);

for idFB = 1:numFloquetPoints_neg
  fprintf('%d sur %d\n', idFB, numFloquetPoints_neg);
  
  FloquetVar = FloquetPoints_neg(idFB);

  % The FE matrix is a linear combination of the elementary pieces
  AAneg = mat_gradu_gradv_neg + 1i * FloquetVar * mat_vecZu_gradv_neg - 1i * FloquetVar * mat_gradu_vecZv_neg + FloquetVar * FloquetVar * mat_vecZu_vecZv_neg  - (opts.omega^2) * mat_u_v_neg;

  % Compute the Floquet-Bloch transform of the solution
  [~, sol_neg_data{idFB}.E0, sol_neg_data{idFB}.E1, sol_neg_data{idFB}.R, sol_neg_data{idFB}.D, BCstruct_neg, TFBlambdaNeg{idFB}] = PeriodicHalfGuideBVP(mesh2Dneg, -1, semiInfiniteDirection, AAneg, BCstruct_neg, numCellsSemiInfinite_neg, opts);
end

%% Positive half-space DtN
TFBphiPosVec = zeros(mesh2Dpos.domain('xmin').numPoints, numBasisInt, numFloquetPoints_pos);
for idI = 1:numBasisInt
  TFBphiPosVec(:, idI, :) = FEPack.tools.BlochTransform(mesh2Dpos.points(mesh2Dpos.domain('xmin').IdPoints, 2),...
                                           FloquetPoints_pos, @(x) int_basis_functions(x, idI),...
                                           1, period_pos, 1000);
end

% TFBphiPos = TFB_Lagrange_P1_1D(mesh2Dint, mesh2Dpos.points(mesh2Dpos.domain('xmin').IdPoints, 2), FloquetPoints_pos, period_pos);
TFBphiPos = cell(numFloquetPoints_pos, 1);
lambda_pos = zeros(numBasisInt);
wpos = (2*pi / period_pos) / (numFloquetPoints_pos - 1);

for idFB = 1:numFloquetPoints_pos
  TFBphiPos{idFB} = BCstruct_pos.spB0.FE_to_spectral * TFBphiPosVec(:, :, idFB);
  TFBlambdaPos{idFB} = BCstruct_pos.spB0.massmat * TFBlambdaPos{idFB};

  lambda_pos = lambda_pos + wpos * TFBphiPos{idFB}' * TFBlambdaPos{idFB} * TFBphiPos{idFB};
end

%% Negative half-space DtN
TFBphiNegVec = zeros(mesh2Dneg.domain('xmin').numPoints, numBasisInt, numFloquetPoints_neg);
for idI = 1:numBasisInt
  TFBphiNegVec(:, idI, :) = FEPack.tools.BlochTransform(mesh2Dneg.points(mesh2Dneg.domain('xmin').IdPoints, 2),...
                                           FloquetPoints_neg, @(x) int_basis_functions(x, idI),...
                                           1, period_neg, 1000);
end

% TFBphiNeg = TFB_Lagrange_P1_1D(mesh2Dint, mesh2Dneg.points(mesh2Dneg.domain('xmin').IdPoints, 2), FloquetPoints_neg, period_neg); 
TFBphiNeg = cell(numFloquetPoints_neg, 1);
lambda_neg = zeros(numBasisInt);
wneg = (2*pi / period_neg) / (numFloquetPoints_neg - 1);

for idFB = 1:numFloquetPoints_neg
  TFBphiNeg{idFB} = BCstruct_neg.spB0.FE_to_spectral * TFBphiNegVec(:, :, idFB);
  TFBlambdaNeg{idFB} = BCstruct_neg.spB0.massmat * TFBlambdaNeg{idFB};

  lambda_neg = lambda_neg + wneg * TFBphiNeg{idFB}' * TFBlambdaNeg{idFB} * TFBphiNeg{idFB};
end

%% Solve the integral equation on the interface
mat_G_v_int = spBint.projmat * Gint(mesh2Dint.points);
% -------------------------------------------------------------------- %
% The minus sign comes from the definition of the Lambda when they are %
% the DtN operators (see PeriodicHalfGuideBVP.m)                       %
% -------------------------------------------------------------------- %
trace_solution = -(lambda_pos + lambda_neg) \ mat_G_v_int;             %
% -------------------------------------------------------------------- %

%% Deduce the FB transform of the solution
% FB transform of the solution's trace with period_pos
TFB_solution_pos = cell(numFloquetPoints_pos, 1);
for idFB = 1:numFloquetPoints_pos

  TFB_solution_pos{idFB} = zeros(mesh2Dpos.numPoints, numCellsSemiInfinite_pos);
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
TFB_solution_neg = cell(numFloquetPoints_neg, 1);
for idFB = 1:numFloquetPoints_neg

  TFB_solution_neg{idFB} = zeros(mesh2Dneg.numPoints, numCellsSemiInfinite_neg);
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

%% Compute the half-space solution cell by cell
% Positive side
numCells = [1 1];
numCells(infiniteDirection) = 2*numCellsInfinite_pos;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_pos;
Nu = prod(numCells);
[I1, I2] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2]; % 2-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite_pos' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period_pos) ./ (numFloquetPoints_pos - 1)); % 1-by-1
U.positive = zeros(mesh2Dpos.numPoints, Nu);

for idFB = 1:numFloquetPoints_pos
  FloquetVar = FloquetPoints_pos(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh2Dpos.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period_pos); % 1-by-Nu
  U_TFB = TFB_solution_pos{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U.positive = U.positive + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U.positive = U.positive * sqrt(period_pos / (2*pi));

% Negative side
numCells = [1 1];
numCells(infiniteDirection) = 2*numCellsInfinite_neg;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_neg;
Nu = prod(numCells);
[I1, I2] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2]; % 2-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite_neg' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period_neg) ./ (numFloquetPoints_neg - 1)); % 1-by-1
U.negative = zeros(mesh2Dneg.numPoints, Nu);

for idFB = 1:numFloquetPoints_neg
  FloquetVar = FloquetPoints_neg(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh2Dneg.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period_neg); % 1-by-Nu
  U_TFB = TFB_solution_neg{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U.negative = U.negative + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U.negative = U.negative * sqrt(period_neg / (2*pi));

%% Plot the solution
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
for idS = 1:numCellsSemiInfinite_pos
  for idI = 1:2*numCellsInfinite_pos
    Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite_pos], idS, idI);
    X = mesh2Dpos.points(:, 1) + (idS - 1);
    Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite_pos - 1) * period_pos;

    trisurf(mesh2Dpos.triangles, X, Y, real(U.positive(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end
%
for idS = 1:numCellsSemiInfinite_neg
  for idI = 1:2*numCellsInfinite_neg
    Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite_neg], idS, idI);
    X = mesh2Dneg.points(:, 1) - (idS - 1);
    Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite_neg - 1) * period_neg;
    trisurf(mesh2Dneg.triangles, X, Y, real(U.negative(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

axis([-4, 4, -4, 4]);
% xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);

