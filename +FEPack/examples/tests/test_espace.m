clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.1i;
theta = pi/3;
opts.cutvec = [cos(theta), sin(theta)];
coSemiInf = 1;
coInf = 2;

mu_pos = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho_pos = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));

mu_neg = @(x)  1 + 0.5*sin(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));% ones(size(x, 1), 1);
rho_neg = @(x) 1 + 0.25*sin(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2)); % ones(size(x, 1), 1);

mu_int = @(x) (x(:, coSemiInf) >= 0) .* mu_pos(x) + (x(:, coSemiInf) < 0) .* mu_neg(x);
rho_int = @(x) (x(:, coSemiInf) >= 0) .* rho_pos(x) + (x(:, coSemiInf) < 0) .* rho_neg(x);

O = [0.5, 0.5]; O(coSemiInf) = -1;
F = @(x) FEPack.tools.cutoff(sqrt((x(:, 1)-O(1)).^2 + (x(:, 2)-O(2)).^2), -0.3, 0.3);

intBounds.pos = +0.5;
intBounds.neg = -2;
structmesh = 0;
basis_functions = 'Lagrange';
u = pdes.PDEObject; v = dual(u);
volBilinearIntg = @(muco, rhoco) (muco * grad(u)) * grad(v) - (opts.omega^2) * ((rhoco*id(u))*id(v));
% volBilinearIntg = @(muco, rhoco) (muco * gradDir(u, opts.cutvec)) * gradDir(v, opts.cutvec) - (opts.omega^2) * ((rhoco*id(u))*id(v));
plot_coefficients = false;

N = 32;

%% Parameters for the positive half-guide
%  //////////////////////////////////////
BB = [0, 1; 0, 1]; BB(coSemiInf, 2) = +1;
mesh_pos = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{2*coSemiInf});
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{2*coSemiInf-1});
else
  FourierIds = [0 0]; FourierIds(3-coSemiInf) = floor(N/4);
  BCstruct_pos.spB0 = spaces.FourierBasis(mesh_pos.domains{2*coSemiInf}, FourierIds);
  BCstruct_pos.spB1 = spaces.FourierBasis(mesh_pos.domains{2*coSemiInf-1}, FourierIds);
end

BCstruct_pos.BCdu = 0.0;
BCstruct_pos.BCu = 1.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_pos.representation = '';
numCells_pos = 4;
volBilinearIntg_pos = volBilinearIntg(@(x)  mu_pos(x + sparse(1:size(x, 1), coSemiInf, intBounds.pos, size(x, 1), size(x, 2))),...
                                      @(x) rho_pos(x + sparse(1:size(x, 1), coSemiInf, intBounds.pos, size(x, 1), size(x, 2))));

%% Parameters for the negative half-guide
%  //////////////////////////////////////
BB = [0, 1; 0, 1]; BB(coSemiInf, 2) = -1;
mesh_neg = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{2*coSemiInf});
  BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{2*coSemiInf-1});
else
  FourierIds = [0 0]; FourierIds(3-coSemiInf) = floor(N/4);
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh_neg.domains{2*coSemiInf}, FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh_neg.domains{2*coSemiInf-1}, FourierIds);
end

BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_neg.representation = '';
numCells_neg = 4;
volBilinearIntg_neg = volBilinearIntg(@(x)  mu_neg(x + sparse(1:size(x, 1), coSemiInf, intBounds.neg, size(x, 1), size(x, 2))),...
                                      @(x) rho_neg(x + sparse(1:size(x, 1), coSemiInf, intBounds.neg, size(x, 1), size(x, 2))));

%% Parameters for the interior problem
%  //////////////////////////////////////
BB = [0, 1; 0, 1]; BB(coSemiInf, 1) = intBounds.neg; BB(coSemiInf, 2) = intBounds.pos;
mesh_int = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);

if strcmpi(basis_functions, 'Lagrange')
  spBint_neg = FEPack.spaces.PeriodicLagrangeBasis(mesh_int.domains{2*coSemiInf});
  spBint_pos = FEPack.spaces.PeriodicLagrangeBasis(mesh_int.domains{2*coSemiInf-1});
else
  FourierIds = [0 0]; FourierIds(3-coSemiInf) = floor(N/4);
  spBint_pos = spaces.FourierBasis(mesh_int.domains{2*coSemiInf}, FourierIds);
  spBint_neg = spaces.FourierBasis(mesh_int.domains{2*coSemiInf-1}, FourierIds);
end

volBilinearIntg_int = volBilinearIntg(mu_int, rho_int);
volLinearIntg = F * id(v);

numCellsSemiInfinite_pos = 4;
numCellsSemiInfinite_neg = 4;
numCellsInfinite = 3;
numFloquetPoints = 100;

%% Plot the coefficients and the source term
if (true)%(plot_coefficients)
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');
  numCells = [1 1];
  numCells(coInf) = 2*numCellsInfinite;
  numCells(coSemiInf) = 2*numCellsSemiInfinite_pos;
  Icells = [1 1];
  orientation = 1;
  for idS = 1:numCells(coSemiInf)
    for idI = 1:numCells(coInf)
      Icells(coInf) = idI;
      Icells(coSemiInf) = idS;

      idcell = sub2ind(numCells, Icells(1), Icells(2));

      X = mesh_pos.points(:, 1) + (coSemiInf == 1) * orientation * (idS - numCellsSemiInfinite_pos - 1) + (coInf == 1) * (idI - numCellsInfinite - 1);
      Y = mesh_pos.points(:, 2) + (coSemiInf == 2) * orientation * (idS - numCellsSemiInfinite_pos - 1) + (coInf == 2) * (idI - numCellsInfinite - 1);
      figure(1);
      trisurf(mesh_pos.triangles, X, Y, mu_int([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

      figure(2);
      trisurf(mesh_pos.triangles, X, Y, rho_int([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

      figure(3);
      trisurf(mesh_pos.triangles, X, Y, F([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end
end

%%
% Compute guide solution
U = PeriodicSpaceBVP(coSemiInf, coInf,...
                     volBilinearIntg_pos, mesh_pos, BCstruct_pos, numCellsSemiInfinite_pos,...
                     volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCellsSemiInfinite_neg,...
                     volBilinearIntg_int, volLinearIntg, mesh_int, spBint_pos, spBint_neg,...
                     numCellsInfinite, numFloquetPoints, opts);

%% Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_pos;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = 1;
for idS = 1:numCells(coSemiInf)
  for idI = 1:numCells(coInf)
    Icells(coInf) = idI;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_pos.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1) + intBounds.pos) + (coInf == 1) * (idI - numCellsInfinite - 1);
    Y = mesh_pos.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1) + intBounds.pos) + (coInf == 2) * (idI - numCellsInfinite - 1);

    trisurf(mesh_pos.triangles, X, Y, real(U.positive(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_neg;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = -1;
for idS = 1:numCells(coSemiInf)
  for idI = 1:numCells(coInf)
    Icells(coInf) = idI;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_neg.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1) + intBounds.neg) + (coInf == 1) * (idI - numCellsInfinite - 1);
    Y = mesh_neg.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1) + intBounds.neg) + (coInf == 2) * (idI - numCellsInfinite - 1);

    trisurf(mesh_neg.triangles, X, Y, real(U.negative(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

numCells = [1 1];
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = 1;
for idS = 1:numCells(coSemiInf)
  for idI = 1:numCells(coInf)
    Icells(coInf) = idI;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_int.points(:, 1) + (coSemiInf == 1) * orientation * (idS - 1) + (coInf == 1) * (idI - numCellsInfinite - 1);
    Y = mesh_int.points(:, 2) + (coSemiInf == 2) * orientation * (idS - 1) + (coInf == 2) * (idI - numCellsInfinite - 1);

    trisurf(mesh_int.triangles, X, Y, real(U.interior(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

BBS = [intBounds.neg - numCellsSemiInfinite_neg + 1, intBounds.pos + numCellsSemiInfinite_pos - 1];
BBI = [-numCellsInfinite, numCellsInfinite];

xlim((coInf == 1) * BBI + (coSemiInf == 1) * BBS);
ylim((coInf == 2) * BBI + (coSemiInf == 2) * BBS);
