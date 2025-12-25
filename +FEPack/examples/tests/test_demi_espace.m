clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.1i;
theta = pi/3;
opts.cutvec = [cos(theta), sin(theta)];
orientation = 1;
coSemiInf = 2;
coInf = 1;
mu = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));

N = 32;
BB = [0, 1; 0, 1]; BB(coSemiInf, 2) = orientation;
mesh = meshes.MeshRectangle(1, BB(1, :), BB(2, :), N, N);
cell = mesh.domain('volumic');
Sigma0 = mesh.domains{2*coSemiInf};
Sigma1 = mesh.domains{2*coSemiInf-1};
%
u = pdes.PDEObject; v = dual(u);
% muGrad = mu;
% volIntg = (muGrad * gradDir(u, opts.cutvec))*gradDir(v, opts.cutvec) - (opts.omega^2) * ((rho*id(u))*id(v));
muGrad = @(x) kron(eye(3), mu(x));
volIntg = (muGrad * grad(u))*grad(v) - (opts.omega^2) * ((rho*id(u))*id(v));

BCstruct.spB0 = FEPack.spaces.PeriodicLagrangeBasis(Sigma0);
BCstruct.spB1 = FEPack.spaces.PeriodicLagrangeBasis(Sigma1);

Nb = BCstruct.spB0.numBasis;
BCstruct.BCdu = 0.0;
BCstruct.BCu = 1.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct.representation = '';
BCstruct.phi = @(x) tools.cutoff(x(:, coInf), -1, 1, 0.5);

numCellsSemiInfinite = 4;
numCellsInfinite = 3;
numFloquetPoints = 100;

U = PeriodicHalfSpaceBVP(mesh, orientation, coSemiInf, coInf, volIntg, ...
                         BCstruct, numCellsSemiInfinite, numCellsInfinite,...
                         numFloquetPoints, opts);
U = full(U);
% profile VIEWER
%
%%  Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite;
numCells(coInf) = 2*numCellsInfinite;

Icells = [1 1];

for idS = 1:numCells(coSemiInf)
  for idI = 1:numCells(coInf)

    Icells(coInf) = idI;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh.points(:, 1) + (coSemiInf == 1) * orientation * (idS - 1) + (coInf == 1) * (idI - numCellsInfinite - 1);
    Y = mesh.points(:, 2) + (coSemiInf == 2) * orientation * (idS - 1) + (coInf == 2) * (idI - numCellsInfinite - 1);

    trisurf(mesh.triangles, X, Y, real(U(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

  end
end

%%
% figure
% phi = @(x) FEPack.tools.cutoff(x, -1, 1, 1);
% tX = linspace(-3, 3, 128);
% tY = linspace(0, 4, 128);
% [X, Y] = meshgrid(tX, tY);
%
% surf(X, Y, phi(X - Y/tan(theta))); shading interp;

%%
% x = linspace(-8, 8, 1024)';
% phi = @(x) FEPack.tools.cutoff(x(:, 1), -1, 1, 1);
% 
% TFBphi = @(x, k) FEPack.tools.BlochTransform(x, k, phi, 1);
% Nk = 100;
% K = linspace(-pi, pi, Nk)';
% phir = @(x) (1/sqrt(2*pi)) * (2*pi/(Nk-1)) * diag(TFBphi(x, K) * exp(+1i*K*x'));
% 
% plot(x, real(phi(x))); hold on;
% plot(x, real(phir(x))); hold on;
