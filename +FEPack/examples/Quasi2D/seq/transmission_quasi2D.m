clear; clc;
warning('Fichiers dans outputs/ supprimés.');
!rm outputs/*
%%
import FEPack.*
% profile OFF
% profile ON

%% Problem-related variables
omega = 8 + 0.5i;
opts.period = 1;
opts.verbose = 0;
problem_setting = 'A'; % 'A' or 'B'

if strcmpi(problem_setting, 'A')

  % 2D coefficients
  period_posFun = 1;
  period_negFun = 0.5*sqrt(2);%*sqrt(2);

  mu2Dpos =  @(x) ones(size(x, 1), 1); % 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2]);
  rho2Dpos = @(x) ones(size(x, 1), 1); % 0.5 + tools.FUN_per_cutoff_cuboid(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1); % 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)/period_negFun);
  rho2Dneg = @(x) ones(size(x, 1), 1); % 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2)/period_negFun);

  % Cut vector
  period_pos = period_posFun;
  period_neg = period_negFun;

  opts.cutvec = [period_pos; period_neg];

  % 3D coefficients
  mu3Dpos = @(x) mu2Dpos([x(:, 1), period_pos * x(:, 2)]);
  rho3Dpos = @(x) rho2Dpos([x(:, 1), period_pos * x(:, 2)]);
  mu3Dneg = @(x) mu2Dneg([x(:, 1), period_neg * x(:, 3)]);
  rho3Dneg = @(x) rho2Dneg([x(:, 1), period_neg * x(:, 3)]);
  
else

  % 2D coefficients
  vecperFun = [-sqrt(2), 1]; % [-sqrt(2), 1];
  
  mu2Dpos = @(x) 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2]);
  rho2Dpos = @(x) 0.5 + tools.FUN_per_cutoff_cuboid(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1);
  rho2Dneg = @(x) ones(size(x, 1), 1);

  % Cut vector
  vecper = vecperFun;
  opts.cutvec = [1/vecper(2); -vecper(1)/vecper(2)];
  
  % 3D coefficients
  fun3D = @(fun2D, x) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
  mu3Dpos = @(x) fun3D(mu2Dpos, x);
  mu3Dneg = @(x) fun3D(mu2Dneg, x);
  rho3Dpos = @(x) fun3D(rho2Dpos, x);
  rho3Dneg = @(x) fun3D(rho2Dneg, x);

end

% Cut matrix and cut slope
opts.cutmat = [[1; 0; 0], [0; opts.cutvec]];

% Jump data
% alpha_G = 3;
% eps_G = 1e-8;
% supp_G = -log(eps_G) / alpha_G;
% G = @(x) exp(-alpha_G * x(:, 2).^2) .* (abs(x(:, 2)) <= supp_G);
G = @(x) FEPack.tools.cutoff(x(:, 2), -0.5, 0.5);
G3D = @(x) G([zeros(size(x, 1), 1), x(:, 2)/opts.cutvec(1), zeros(size(x, 1), 1)]);

% Numbers of cells
numCellsSemiInfinite_pos = 6;
numCellsSemiInfinite_neg = 6;
opts.numCellsInfinite = 6;
opts.numFloquetPoints = 50;
opts.FloquetPoints = linspace(-pi/opts.period, pi/opts.period, opts.numFloquetPoints);

%% Meshes
struct_mesh = 0;

numNodesX = 8;
numNodesY = 8;
numNodesZ = 8;
ratio_s = 4; % more nodes for the transverse direction

meshXYpos = meshes.MeshRectangle(struct_mesh, [0  1], [0 1/opts.cutvec(1)], numNodesX, ceil(numNodesY/opts.cutvec(1)));
meshXYneg = meshes.MeshRectangle(struct_mesh, [0 -1], [0 1/opts.cutvec(1)], numNodesX, ceil(numNodesY/opts.cutvec(1)));
meshZXpos = meshes.MeshRectangle(struct_mesh, [0 1], [0  1], numNodesZ, numNodesX);
meshZXneg = meshes.MeshRectangle(struct_mesh, [0 1], [0 -1], numNodesZ, numNodesX);
meshYZ = meshes.MeshRectangle(struct_mesh, [0 1], [0 1], numNodesY, ratio_s*numNodesZ);

%% Basis functions and boundary conditions
% Positive side
basis_functions = 'Fourier';
if strcmpi(basis_functions, 'Lagrange')
  BCstruct.pos.spBX = FEPack.spaces.PeriodicLagrangeBasis(meshYZ.domain('volumic')); % X = cst
else
  BCstruct.pos.spBX = spaces.FourierBasis(meshYZ.domain('volumic'), [floor(numNodesY/4), floor(numNodesZ/4)]);
end
% Particular Lagrange for basis functons on Y = cst
ecs_fun = @(u, dom) ((((u|dom.mesh.domain('xmin')) - (u|dom.mesh.domain('xmax'))) == 0.0) &...      % Periodic on Z
                      ((u|dom.mesh.domain('ymin')) == 0.0) & ((u|dom.mesh.domain('ymax')) == 0.0)); % Dirichlet on X

BCstruct.pos.spBY = FEPack.spaces.EssentialBCLagrangeBasis(meshZXpos.domain('volumic'), ecs_fun); % Y = cst

% Other parameters
BCstruct.pos.BCdu = 0.0;
BCstruct.pos.BCu = 1.0;
BCstruct.pos.representation = 'projection';

% Negative side
BCstruct.neg = BCstruct.pos;

% Particular Lagrange for basis functons on Y = cst
ecs_fun = @(u, dom) ((((u|dom.mesh.domain('xmin')) - (u|dom.mesh.domain('xmax'))) == 0.0) &...      % Periodic on Z
                      ((u|dom.mesh.domain('ymin')) == 0.0) & ((u|dom.mesh.domain('ymax')) == 0.0)); % Dirichlet on X

BCstruct.neg.spBY = FEPack.spaces.EssentialBCLagrangeBasis(meshZXneg.domain('volumic'), ecs_fun); % Y = cst

%% Compute DtN operators
opts.suffix = 'pos';
quasi2DHalfGuideDirichlet(+1, meshXYpos, meshYZ, omega, mu3Dpos, rho3Dpos, BCstruct.pos, opts);

%%
opts.suffix = 'neg';
quasi2DHalfGuideDirichlet(-1, meshXYneg, meshYZ, omega, mu3Dneg, rho3Dneg, BCstruct.neg, opts);

%% Interface equation
fprintf('<strong>Equation d''interface.</strong>\n');
for idFB = 1:opts.numFloquetPoints
  fprintf('%d sur %d\n', idFB, opts.numFloquetPoints);
  FloquetVar = opts.FloquetPoints(idFB);

  % DtN operators
  solguidePos = load(['outputs/half_guide_sol_pos_Floquet_', num2str(idFB)]);
  solguideNeg = load(['outputs/half_guide_sol_neg_Floquet_', num2str(idFB)]);
  
  Lambda_pos = solguidePos.Lambda;
  Lambda_neg = solguideNeg.Lambda;

  % Rhs
  jumpData_FB = @(x) FEPack.tools.BlochTransform(x, FloquetVar, G3D, 2, opts.period);
  GG = BCstruct.pos.spBX.FE_to_spectral * jumpData_FB(meshYZ.points);
    
  % The minus sign comes from the definition of the Lambda
  soltrace = -(Lambda_pos + Lambda_neg) \ GG;
  
  % Save the trace of the solution
  save(['outputs/sol_trace_Floquet_', num2str(idFB)], 'soltrace', '-v7.3');
end

%% Construct solution
opts.suffix = 'pos';
opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
construct_solution(meshXYpos, meshYZ, BCstruct.pos, opts);

opts.suffix = 'neg';
opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
construct_solution(meshXYneg, meshYZ, BCstruct.neg, opts);

% Plot solution
mafig = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idY = 1:2*opts.numCellsInfinite
  sol = load(['outputs/solution_pos_Y_', num2str(idY)]);

  for idX = 1:numCellsSemiInfinite_pos
    X = meshXYpos.points(:, 1) + (idX - 1);
    Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) / opts.cutvec(1);

    trisurf(meshXYpos.triangles, X, Y, real(sol.Usol(:, idX)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end
%
%
for idY = 1:2*opts.numCellsInfinite
  sol = load(['outputs/solution_neg_Y_', num2str(idY)]);

  for idX = 1:numCellsSemiInfinite_neg
    X = meshXYneg.points(:, 1) - (idX - 1);
    Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) / opts.cutvec(1);

    trisurf(meshXYneg.triangles, X, Y, real(sol.Usol(:, idX)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
savefig(mafig, 'outputs/solution', 'compact');
 