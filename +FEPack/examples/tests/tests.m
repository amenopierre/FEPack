% clear; clc;
% import FEPack.*
% 
% % profile ON
% opts.omega = 8 + 0.1i;
% theta = pi/3;
% opts.cutvec = [cos(theta), sin(theta)];
% orientation = 1;
% coSemiInf = 2;
% coInf = [1 3];
% mu = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2)) + 0.25*cos(2*pi*x(:, 3));
% rho = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2)).*sin(2*pi*x(:, 3));
% 
% N = 8;
% BB = [0, 1; 0, 1; 0, 1]; BB(coSemiInf, 2) = orientation;
% mesh = meshes.MeshCuboid(1, BB(1, :), BB(2, :), BB(3, :), N, N);
% cellule = mesh.domain('volumic');
% Sigma0 = mesh.domains{2*coSemiInf};
% Sigma1 = mesh.domains{2*coSemiInf-1};
% %
% u = pdes.PDEObject; v = dual(u);
% % muGrad = mu;
% % volIntg = (muGrad * gradDir(u, opts.cutvec))*gradDir(v, opts.cutvec) - (opts.omega^2) * ((rho*id(u))*id(v));
% muGrad = @(x) kron(eye(3), mu(x));
% volIntg = (muGrad * grad(u))*grad(v) - (opts.omega^2) * ((rho*id(u))*id(v));
% 
% BCstruct.spB0 = FEPack.spaces.PeriodicLagrangeBasis(Sigma0);
% BCstruct.spB1 = FEPack.spaces.PeriodicLagrangeBasis(Sigma1);
% 
% Nb = BCstruct.spB0.numBasis;
% BCstruct.BCdu = 0.0;
% BCstruct.BCu = 1.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
% BCstruct.representation = '';
% BCstruct.phi = @(x) tools.cutoff(x(:, 1), -1, 1, 0.5);
% 
% numCellsSemiInfinite = 4;
% numCellsInfinite = [3 3];
% numFloquetPoints = [5 5];
% 
% U = PeriodicHalfSpaceBVP(mesh, orientation, coSemiInf, coInf, volIntg, ...
%                          BCstruct, numCellsSemiInfinite, numCellsInfinite,...
%                          numFloquetPoints, opts);
% U = full(U);

%%
% N2D = 32;
% mesh2D = meshes.MeshCuboid(1, BB(1, :), BB(2, :), BB(3, :), N, N);

% profile VIEWER
%
% %%  Plot U
% figure;
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
%
% numCells = [1 1];
% numCells(coSemiInf) = numCellsSemiInfinite;
% numCells(coInf) = 2*numCellsInfinite;
%
% Icells = [1 1];
%
% for idS = 1:numCells(coSemiInf)
%   for idI = 1:numCells(coInf)
%
%     Icells(coInf) = idI;
%     Icells(coSemiInf) = idS;
%
%     idcell = sub2ind(numCells, Icells(1), Icells(2));
%
%     X = mesh.points(:, 1) + (coSemiInf == 1) * orientation * (idS - 1) + (coInf == 1) * (idI - numCellsInfinite - 1);
%     Y = mesh.points(:, 2) + (coSemiInf == 2) * orientation * (idS - 1) + (coInf == 2) * (idI - numCellsInfinite - 1);
%
%     trisurf(mesh.triangles, X, Y, real(U(:, idcell)));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%
%   end
% end

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

%%
% x = linspace(-3, 3, 1024)';
% h = 0.5;
% phi = @(x) (-abs(x(:, 1)/h)+1).*(abs(x(:, 1)) < h);
% TFBphi = @(x, k) FEPack.tools.BlochTransform(x, k, @(x) phi(x), 1, 1, 100, 'periodic');
% 
% % figure;
% % plot(x, phi(x));
% 
% figure;
% Nk = 1024;
% K = linspace(-pi, pi, Nk)';
% TFBphiVec = TFBphi(x, K);
% pause;
% for i = 1:Nk
%   subplot(2, 1, 1);
%   plot(x, real(TFBphiVec(:, i)));
%   ylim([-0.5, 0.5]);
%   title(['k = ', num2str(K(i), '%d')], 'interpreter', 'latex');
% 
%   subplot(2, 1, 2);
%   plot(x, imag(TFBphiVec(:, i)));
%   ylim([-0.5, 0.5]);
% 
%   pause(0.1);
%   hold off;
% end
% 
% %%
% clear; clc;
% import FEPack.*
% 
% mesh = FEPack.meshes.MeshRectangle(0, [0, 1], [0, 1], 32, 32);
% Omega = mesh.domain('volumic');
% Sigma = mesh.domain('xmin');
% 
% u = FEPack.pdes.PDEObject;
% v = dual(u);
% 
% AA = FEPack.pdes.Form.intg(Omega, (grad(u))*grad(v) + id(u)*id(v));


%%
clear; clc;
import FEPack.*
N = 32;

omega = @(x) ones(size(x, 1), 1);
% mu = @(x) 1.5 + cos(2*pi*x(:, 1));
rho = @(x) ones(size(x, 1), 1);
% rho = @(x) 1 + 0.5*sin(2*pi*x(:, 1));

mesh = FEPack.meshes.MeshSegment('uniform', 0, 1, N);
Omega = mesh.domain('volumic');
x0 = mesh.domain('xmin');
x1 = mesh.domain('xmax');

u = FEPack.pdes.PDEObject;
v = dual(u);

MM = FEPack.pdes.Form.intg(Omega, rho * u * v);
KK = FEPack.pdes.Form.intg(Omega, (omega * grad(u))*grad(v));

ecsDir = ((u|x0) == 0) & ((u|x1) == 0); ecsDir.applyEcs;
ecsPer = ((u|x0) - (u|x1) == 0); ecsPer.applyEcs;

[phiDir, lambdaDir] = eigs(ecsDir.P * KK * ecsDir.P', ecsDir.P * MM * ecsDir.P', N-2, 'sm');
[phiPer, lambdaPer] = eigs(ecsPer.P * KK * ecsPer.P', ecsPer.P * MM * ecsPer.P', N-2, 'sm');
lambdaDir = diag(lambdaDir);
lambdaPer = diag(lambdaPer);
phiDir = ecsDir.P' * phiDir;
phiPer = ecsPer.P' * phiPer;

plot((lambdaDir), 'bx', 'markersize', 10); hold on;
plot((lambdaPer), 'ro', 'markersize', 10); hold on;
% plot((lambdaDir), zeros(size(lambdaDir)), 'bx', 'markersize', 10); hold on;
% plot((lambdaPer), zeros(size(lambdaPer)), 'ro', 'markersize', 10); hold on;
%%
omega = linspace(0, 10, 1024);

Nev = 1000;
L = 1;
k = (1:Nev)';
lambda = (k*pi/L).^2;
alpha0 = (k*pi/L)*sqrt(2/L);
alpha1 = -alpha0.*((-1).^k);
zeta = 1 ./ (omega.^2 - lambda);

t00 = (alpha0.*alpha0)' * zeta;
t10 = (alpha1.*alpha0)' * zeta;
t01 = (alpha0.*alpha1)' * zeta;
t11 = (alpha1.*alpha1)' * zeta;


plot(omega, t11);
hold on;
plot(omega, omega./tan(omega*L));
% plot(omega, -omega./sin(omega*L))
hold off;

%%
% mesh = FEPack.meshes.MeshSegment('uniform', 0, 1, 32);
% A = eye(2);
% D = zeros(2);
% b = ones(2);
% 
% sol = FEPack.solver.CellBVP;
% sol.initialize(2, 'numEdgeNodes', [8 8]);

%%
% findArea(3, 4)
% 
% function a = findArea(hola, width,varargin)
%    defaultHeight = 1;
%    defaultUnits = 'inches';
%    defaultShape = 'rectangle';
%    expectedShapes = {'square','rectangle','parallelogram'};
% 
%    p = inputParser;
%    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%    addRequired(p,'width', validScalarPosNum);
%    addOptional(p,'height',defaultHeight,validScalarPosNum);
% %    addParameter(p,'units',defaultUnits,@isstring);
% %    addParameter(p,'shape',defaultShape,...
% %                  @(x) any(validatestring(x,expectedShapes)));
%    parse(p, hola, width,varargin{:});
%    p
% %    
% %    a = p.Results.width*p.Results.height; 
% end
