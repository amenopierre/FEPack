clear; clc;
%%
import FEPack.*
% profile OFF
% profile ON

%% Problem-related variables
omega = 8 + 0.25i;
opts.omega = omega;
period = 1;
opts.verbose = 0;
problem_setting = 'B'; % 'A' or 'B'

if strcmpi(problem_setting, 'A')

  % 2D coefficients
  period_posFun = 1;
  period_negFun = sqrt(2);

  mu2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2]);
  rho2Dpos = @(x) 0.5 + perCutoffCuboid(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu2Dneg  = @(x) 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)/period_negFun);
  rho2Dneg = @(x) 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2)/period_negFun);

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
  
  mu2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2]);
  rho2Dpos = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
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
G = @(x) FEPack.tools.cutoff(x(:, 2), -0.5, 0.5);
G3D = @(x) G([zeros(size(x, 1), 1), x(:, 2)/cutvec(1), zeros(size(x, 1), 1)]);

% (semi-)infinite directions and numbers of cells
semiInfiniteDirection = 1;
infiniteDirection = 2;
numCellsSemiInfinite_pos = 8;
numCellsSemiInfinite_neg = 8;
numCellsInfinite = 8;
numFloquetPoints = 64;

%% Mesh
pregenerate_mesh = 1;
struct_mesh = 1;
numNodes2D = 60;
numNodes3D = 40;

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
  
  error;
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
fprintf('Computing mat_gradu_gradv_pos...\n');
mat_gradu_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_gradv(mu3Dpos));
fprintf('Computing mat_gradu_gradv_pos...\n');
mat_gradu_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_vec1v(mu3Dpos));
fprintf('Computing mat_gradu_gradv_pos...\n');
mat_vec1u_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_gradv(mu3Dpos));
mat_vec1u_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_vec1v(mu3Dpos));
mat_u_v_pos         = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'),        u_v(rho3Dpos));

% Negative side
mat_gradu_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_gradv(mu3Dneg));
mat_gradu_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_vec1v(mu3Dneg));
mat_vec1u_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_gradv(mu3Dneg));
mat_vec1u_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_vec1v(mu3Dneg));
mat_u_v_neg         = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'),        u_v(rho3Dneg));

%% Boundary conditions
basis_functions = 'Lagrange';
if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmin'));
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmax'));
else
  FourierIds = [0, numNodes3D/2, 0];
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
  FourierIds = [0, numNodes3D/2, 0];
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh3Dneg.domain('xmin'), FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh3Dneg.domain('xmax'), FourierIds);
end
BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;
BCstruct_neg.representation = 'projection';

%% Floquet-Bloch transform of the solution
FloquetPoints = linspace(-pi/period, pi/period, numFloquetPoints);
TFBU = cell(numFloquetPoints, 1);

for idFB = 1:numFloquetPoints
  fprintf('%d sur %d\n', idFB, numFloquetPoints);
  
  FloquetVar = FloquetPoints(idFB);

  % The FE matrix is a linear combination of the elementary pieces
  AApos = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vec1u_gradv_pos - 1i * FloquetVar * mat_gradu_vec1v_pos + FloquetVar * FloquetVar * mat_vec1u_vec1v_pos - (omega^2) * mat_u_v_pos;
  AAneg = mat_gradu_gradv_neg + 1i * FloquetVar * mat_vec1u_gradv_neg - 1i * FloquetVar * mat_gradu_vec1v_neg + FloquetVar * FloquetVar * mat_vec1u_vec1v_neg - (omega^2) * mat_u_v_neg;

  % The Floquet-Bloch transform of the boundary data
  jumpData_FB = @(x) FEPack.tools.BlochTransform(x, FloquetVar, G3D, infiniteDirection);

  % Compute the Floquet-Bloch transform of the solution
  TFBU{idFB} = PeriodicGuideJumpBVP(semiInfiniteDirection,...
                                    AApos, mesh3Dpos, BCstruct_pos, numCellsSemiInfinite_pos,...
                                    AAneg, mesh3Dneg, BCstruct_neg, numCellsSemiInfinite_neg,...
                                    jumpData_FB, opts);
end

%% Take the inverse Floquet transform: positive side
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
  FloquetVar = FloquetPoints(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dpos.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFBU{idFB}.positive(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.positive = U3D.positive + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.positive = U3D.positive * sqrt(period / (2*pi));

%% Take the inverse Floquet transform: negative side
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
  FloquetVar = FloquetPoints(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dneg.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFBU{idFB}.negative(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.negative = U3D.negative + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.negative = U3D.negative * sqrt(period / (2*pi));

%% Take the trace: Positive side
N2Dpos = mesh2Dpos.numPoints;
U2D.positive = zeros(N2Dpos, size(U3D.positive, 2));
dom = mesh3Dpos.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_pos*(idI-1)+1):(numCellsSemiInfinite_pos*idI);
  X = mesh2Dpos.points(:, 1);
  Y = mesh2Dpos.points(:, 2); % ones(mesh2Dpos.numPoints, 1);
  Z = FEPack.tools.mymod(cutslope * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.positive(:, IcellY) = reshape(sum(reshape(coos .* U3D.positive(elts, IcellY), dom.dimension+1, []), 1), N2Dpos, []);
end

%% Take the trace: Negative side
N2Dneg = mesh2Dneg.numPoints;
U2D.negative = zeros(N2Dneg, size(U3D.negative, 2));
dom = mesh3Dneg.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_neg*(idI-1)+1):(numCellsSemiInfinite_neg*idI);
  X = mesh2Dneg.points(:, 1);
  Y = mesh2Dneg.points(:, 2); % ones(mesh2Dneg.numPoints, 1);
  Z = FEPack.tools.mymod(cutslope * (Y + idI - numCellsInfinite - 1));
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
figure;%(1);
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
