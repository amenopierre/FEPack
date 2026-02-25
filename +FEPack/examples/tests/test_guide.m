clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.5i;
theta = pi/3;
opts.cutvec = [cos(theta), sin(theta)];
coInf = 1;
mu = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));

O = [0.5, 0.5]; O(coInf) = 0;
F = @(x) FEPack.tools.cutoff(sqrt((x(:, 1)-O(1)).^2 + (x(:, 2)-O(2)).^2), -0.3, 0.3);
intBounds.pos = +1;
intBounds.neg = -1;
structmesh = 0;
basis_functions = 'Fourier';

u = pdes.PDEObject; v = dual(u);
% volIntg = (mu * gradDir(u, opts.cutvec))*gradDir(v, opts.cutvec) - (opts.omega^2) * ((rho*id(u))*id(v));
% volIntg = (muGrad * grad(u))*grad(v) - (opts.omega^2) * ((rho*id(u))*id(v));

volBilinearIntg = @(muco, rhoco) (muco * grad2(u)) * grad2(v) - (opts.omega^2) * ((rhoco*id(u))*id(v));
% volBilinearIntg = @(muco, rhoco) (muco * gradDir(u, opts.cutvec)) * gradDir(v, opts.cutvec) - (opts.omega^2) * ((rhoco*id(u))*id(v));

N = 32;

% Parameters for the positive half-guide
BB = [0, 1; 0, 1]; BB(coInf, 2) = +1;
mesh_pos = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{2*coInf});
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{2*coInf-1});
else
  FourierIds = [0 0]; FourierIds(3-coInf) = floor(N/4);
  BCstruct_pos.spB0 = spaces.FourierBasis(mesh_pos.domains{2*coInf}, FourierIds);
  BCstruct_pos.spB1 = spaces.FourierBasis(mesh_pos.domains{2*coInf-1}, FourierIds);
end

BCstruct_pos.BCdu = 0.0;
BCstruct_pos.BCu = 1.0;% -1i*opts.omega;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_pos.representation = '';
numCells_pos = 4;
volBilinearIntg_pos = volBilinearIntg(@(x)  mu(x + sparse(1:size(x, 1), coInf, intBounds.pos, size(x, 1), size(x, 2))),...
                                      @(x) rho(x + sparse(1:size(x, 1), coInf, intBounds.pos, size(x, 1), size(x, 2))));

% Parameters for the negative half-guide
BB = [0, 1; 0, 1]; BB(coInf, 2) = -1;
mesh_neg = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{2*coInf});
  BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{2*coInf-1});
else
  FourierIds = [0 0]; FourierIds(3-coInf) = floor(N/4);
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh_neg.domains{2*coInf}, FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh_neg.domains{2*coInf-1}, FourierIds);
end

BCstruct_neg.BCdu = 0.0; % 1.0;
BCstruct_neg.BCu = 1.0; % -1i*opts.omega; % @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_neg.representation = '';
numCells_neg = 4;
volBilinearIntg_neg = volBilinearIntg(@(x)  mu(x + sparse(1:size(x, 1), coInf, intBounds.neg, size(x, 1), size(x, 2))),...
                                      @(x) rho(x + sparse(1:size(x, 1), coInf, intBounds.neg, size(x, 1), size(x, 2))));

% Parameters for the interior problem
BB = [0, 1; 0, 1]; BB(coInf, 1) = intBounds.neg; BB(coInf, 2) = intBounds.pos;
mesh_int = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);

if strcmpi(basis_functions, 'Lagrange')
  spBint_neg = FEPack.spaces.PeriodicLagrangeBasis(mesh_int.domains{2*coInf});
  spBint_pos = FEPack.spaces.PeriodicLagrangeBasis(mesh_int.domains{2*coInf-1});
else
  FourierIds = [0 0]; FourierIds(3-coInf) = floor(N/4);
  spBint_neg = spaces.FourierBasis(mesh_int.domains{2*coInf}, FourierIds);
  spBint_pos = spaces.FourierBasis(mesh_int.domains{2*coInf-1}, FourierIds);
end

volBilinearIntg_int = volBilinearIntg(mu, rho);
volLinearIntg = F * id(v);

%%
% Compute guide solution
U = PeriodicGuideBVP(coInf,...
                     volBilinearIntg_pos, mesh_pos, BCstruct_pos, numCells_pos,...
                     volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCells_neg,...
                     volBilinearIntg_int, volLinearIntg, mesh_int, spBint_pos, spBint_neg, opts);

%% Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idI = 1:numCells_pos
  subplot(-coInf+3, coInf, 1);
  trisurf(mesh_pos.triangles, mesh_pos.points(:, 1) + (coInf == 1) * (idI-1+intBounds.pos),...
                              mesh_pos.points(:, 2) + (coInf == 2) * (idI-1+intBounds.pos), real(U.positive(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  xlim([-5, 5]);

  subplot(-coInf+3, coInf, 2);
  trisurf(mesh_pos.triangles, mesh_pos.points(:, 1) + (coInf == 1) * (idI-1+intBounds.pos),...
                              mesh_pos.points(:, 2) + (coInf == 2) * (idI-1+intBounds.pos), imag(U.positive(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  xlim([-5, 5]);
end

for idI = 1:numCells_neg
  subplot(-coInf+3, coInf, 1);
  trisurf(mesh_neg.triangles, mesh_neg.points(:, 1) + (coInf == 1) * (-(idI-1)+intBounds.neg),...
                              mesh_neg.points(:, 2) + (coInf == 2) * (-(idI-1)+intBounds.neg), real(U.negative(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  xlim([-5, 5]);

  subplot(-coInf+3, coInf, 2);
  trisurf(mesh_neg.triangles, mesh_neg.points(:, 1) + (coInf == 1) * (-(idI-1)+intBounds.neg),...
                              mesh_neg.points(:, 2) + (coInf == 2) * (-(idI-1)+intBounds.neg), imag(U.negative(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  xlim([-5, 5]);
end

subplot(-coInf+3, coInf, 1);
trisurf(mesh_int.triangles, mesh_int.points(:, 1), mesh_int.points(:, 2), real(U.interior));
hold on;
view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
xlim([-5, 5]);

subplot(-coInf+3, coInf, 2);
trisurf(mesh_int.triangles, mesh_int.points(:, 1), mesh_int.points(:, 2), imag(U.interior));
hold on;
view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
xlim([-5, 5]);
