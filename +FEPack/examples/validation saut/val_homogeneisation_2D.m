clear; clc;
% homogeneisation_2D(5 + 0.5i, 1,    1);
% homogeneisation_2D(5 + 0.5i, 0.5,  2);
% homogeneisation_2D(5 + 0.5i, 0.25, 3);

homogeneisation_2D(5 + 0.5i);
% homogeneisation_2D(0.25 + 0.25i);
% homogeneisation_2D(0.125 + 0.125i);


function homogeneisation_2D(freq, period, suffix)

  if (nargin < 2)
    period = 1;
  end

  import FEPack.*
  % profile ON
  opts.omega = freq; %0.125 + 0.125i;
  coSemiInf = 1;
  coInf = 2;

  % mu_pos = @(x) 4*ones(size(x, 1), 1);
  % mu_pos = @(x) [ones(size(x, 1), 1), zeros(size(x, 1), 1); zeros(size(x, 1), 1), 2 * ones(size(x, 1), 1)];
  % mu_pos = @(x) [2 + sin(2*pi*x(:, 1)/period), zeros(size(x, 1), 1); zeros(size(x, 1), 1), 4 * ones(size(x, 1), 1)];
  % mu_pos = @(x) (2 + sin(2*pi*x(:, 1))) .* (4 + sin(2*pi*x(:, 2)));
  
  mu_pos = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
  rho_pos = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));

  mu_pos_eff = @(x) 4*ones(size(x, 1), 1);
  % mu_pos_eff = @(x) [ones(size(x, 1), 1), zeros(size(x, 1), 1); zeros(size(x, 1), 1), 2 * ones(size(x, 1), 1)];
  % mu_pos_eff = @(x) [sqrt(3) * ones(size(x, 1), 1), zeros(size(x, 1), 1); zeros(size(x, 1), 1), 4 * ones(size(x, 1), 1)];
  % mu_pos_eff = @(x) [4*sqrt(3) * ones(size(x, 1), 1), zeros(size(x, 1), 1); zeros(size(x, 1), 1), 2*sqrt(15) * ones(size(x, 1), 1)];
  rho_pos_eff = @(x) 2*ones(size(x, 1), 1);

  mu_neg = @(x)  ones(size(x, 1), 1); % 1 + 0.5*sin(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
  rho_neg = @(x) ones(size(x, 1), 1); % 1 + 0.25*sin(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2)); 

  % mu_glo = @(x) (x(:, coSemiInf) >= 0) .* mu_pos(x) + (x(:, coSemiInf) < 0) .* mu_neg(x);
  % rho_glo = @(x) (x(:, coSemiInf) >= 0) .* rho_pos(x) + (x(:, coSemiInf) < 0) .* rho_neg(x);


  G = @(x) FEPack.tools.cutoff(x(:, 2), -1.0, 1.0);

  structmesh = 0;
  basis_functions = 'Lagrange';
  u = pdes.PDEObject; v = dual(u);
  volBilinearIntg = @(muco, rhoco) (muco * grad2(u)) * grad2(v) - (opts.omega^2) * ((rhoco*id(u))*id(v));
  % plot_coefficients = false;

  N = 50;

  %% Parameters for the positive half-guide
  %  //////////////////////////////////////
  % BB = [0, period; 0, period]; BB(coSemiInf, 2) = +period;
  % mesh_pos = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);
  m = load(['pregenMeshes/2D/unstruct_mesh_2D_', int2str(N), '_positive.mat']);
  mesh_pos = m.mesh;

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
  volBilinearIntg_pos = volBilinearIntg(mu_pos, rho_pos);
  volBilinearIntg_pos_eff = volBilinearIntg(mu_pos_eff, rho_pos_eff);

  %% Parameters for the negative half-guide
  %  //////////////////////////////////////
  % BB = [0, period; 0, period]; BB(coSemiInf, 2) = -period;
  % mesh_neg = meshes.MeshRectangle(structmesh, BB(1, :), BB(2, :), N, N);
  m = load(['pregenMeshes/2D/unstruct_mesh_2D_', int2str(N), '_negative_X.mat']);
  mesh_neg = m.mesh;

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
  volBilinearIntg_neg = volBilinearIntg(mu_neg, rho_neg);

  %% Parameters for the interface problem
  %  //////////////////////////////////////
  numCellsSemiInfinite_pos = 5;
  numCellsSemiInfinite_neg = 5;
  numCellsInfinite = 5;
  numFloquetPoints = 50;

  %%
  % Compute guide solution
  U = PeriodicSpaceJumpBVP(coSemiInf, coInf, period,...
                          volBilinearIntg_pos, mesh_pos, BCstruct_pos, numCellsSemiInfinite_pos,...
                          volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCellsSemiInfinite_neg,...
                          G, numCellsInfinite, numFloquetPoints, opts);
  % %
  % %
  % Ueff = PeriodicSpaceJumpBVP(coSemiInf, coInf, period,...
  %                         volBilinearIntg_pos_eff, mesh_pos, BCstruct_pos, numCellsSemiInfinite_pos,...
  %                         volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCellsSemiInfinite_neg,...
  %                         G, numCellsInfinite, numFloquetPoints, opts);

  % %%
  % % Error
  % E.positive = U.positive - Ueff.positive;
  % E.negative = U.negative - Ueff.negative;

  if (nargin >= 3)
    save(['output_', suffix, '.mat']);
  end


  
  %% Plot the error
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

      X = mesh_pos.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1)) * period + (coInf == 1) * (idI - numCellsInfinite - 1) * period;
      Y = mesh_pos.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1)) * period + (coInf == 2) * (idI - numCellsInfinite - 1) * period;

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

      X = mesh_neg.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1)) * period + (coInf == 1) * (idI - numCellsInfinite - 1) * period;
      Y = mesh_neg.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1)) * period + (coInf == 2) * (idI - numCellsInfinite - 1) * period;

      trisurf(mesh_neg.triangles, X, Y, real(U.negative(:, idcell)));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end


  BBS = [-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1];
  BBI = [-numCellsInfinite, numCellsInfinite];

  xlim((coInf == 1) * BBI + (coSemiInf == 1) * BBS);
  ylim((coInf == 2) * BBI + (coSemiInf == 2) * BBS);
  % caxis([-1e-3, 1e-3]);
end