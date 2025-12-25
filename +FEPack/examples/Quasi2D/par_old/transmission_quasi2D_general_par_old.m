function transmission_quasi2D_general_par(idtest)
% clear; clc;
%%
import FEPack.*
profile OFF
% profile ON

%% Problem-related variables
opts.period = 1;
Lint = 10;
opts.verbose = 0;
opts.nomdossier = 'outputs/'; % /!\ Pas oublier le / à la fin
% opts.nomdossier = '/UMA/tmp/amenoagbadji/solution_B_irrational_per_2_G_1/'; % /!\ Pas oublier le / à la fin
nomdossier = opts.nomdossier;

warning(['Fichiers dans ', nomdossier, ' supprimés.']);
system(['rm ', nomdossier, '*']);

mu2Dpos =  @(x) ones(size(x, 1), 1);
mu2Dneg  = @(x) ones(size(x, 1), 1);

if (idtest == 1)

  opts.omega = 8 + 0.25i;
  
  period_posFun = 1;
  period_negFun = sqrt(2);

  rho2Dpos = @(x) 0.5 + perCutoffCircle(x/opts.period, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.4, 0.4]);
  rho2Dneg = @(x) 0.5 + perCutoffCuboid(x/opts.period, [1; 0], [0; period_negFun], [0.5, 0.5], [-0.25, 0.25], [-0.25, 0.25], 0.5, 1);

  vecperfunpos = [0, period_posFun];
  vecperfunneg = [0, period_negFun];

elseif (idtest == 2)

  opts.omega = 8 + 0.25i;
  vecperFun = [cos(3*pi/5), sin(3*pi/5)];
  rho2Dpos = @(x) 0.5 + perCutoffCircle(x/opts.period, [1; 0], vecperFun, [0.5, 0.5], [-0.4, 0.4]);
  rho2Dneg = @(x) ones(size(x, 1), 1);

  vecperfunpos = vecperFun;
  vecperfunneg = [0, 1];

elseif (idtest == 3)

  opts.omega = 8 + 0.25i;
  vecperfunpos = [cos(3*pi/5), sin(3*pi/5)];
  vecperfunneg = [cos(pi/3), sin(pi/3)];

  mu2Dpos =  @(x) ones(size(x, 1), 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1);

  rho2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.4, 0.4]);
  rho2Dneg = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperfunneg, [0.5, 0.5], [-0.25, 0.25], [-0.25, 0.25], 0.5, 1);

elseif (idtest == 4)

  opts.omega = 20 + 0.25i;
  vecperfunpos = [cos(3*pi/5), sin(3*pi/5)];
  vecperfunneg = [cos(pi/3), sin(pi/3)];

  mu2Dpos =  @(x) ones(size(x, 1), 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1);

  rho2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.4, 0.4]);
  rho2Dneg = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperfunneg, [0.5, 0.5], [-0.25, 0.25], [-0.25, 0.25], 0.5, 1);

elseif (idtest == 5)

  opts.omega = 20 + 0.05i;
  vecperfunpos = [cos(3*pi/5), sin(3*pi/5)];
  vecperfunneg = [cos(pi/3), sin(pi/3)];

  mu2Dpos =  @(x) ones(size(x, 1), 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1);

  rho2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.4, 0.4]);
  rho2Dneg = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperfunneg, [0.5, 0.5], [-0.25, 0.25], [-0.25, 0.25], 0.5, 1);

end

% Cut vector
vecperpos = [-vecperfunpos(1), vecperfunpos(2)];
vecperneg = [-vecperfunneg(1), vecperfunneg(2)];

cutvecpos = [1/vecperpos(2); -vecperpos(1)/vecperpos(2)];
cutvecneg = [1/vecperneg(2); -vecperneg(1)/vecperneg(2)];

% Cut matrix and cut slope
cutmat_pos = [[1; 0; 0], [0; cutvecpos]];
cutmat_neg = [[1; 0; 0], [0; cutvecneg]];

cutslopePos = cutvecpos(2) / cutvecpos(1);
cutslopeNeg = cutvecneg(2) / cutvecneg(1);

% 3D coefficients
fun3D = @(fun2D, x, vecper) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
mu3Dpos  = @(x) fun3D(mu2Dpos,  x, vecperpos);
rho3Dpos = @(x) fun3D(rho2Dpos, x, vecperpos);
mu3Dneg  = @(x) fun3D(mu2Dneg,  x, vecperneg);
rho3Dneg = @(x) fun3D(rho2Dneg, x, vecperneg);

% Jump data
Gint = @(x) FEPack.tools.cutoff(x(:, 1), -0.5, 0.5);
Gtrs = @(s) ones(size(s, 1), 1);

% Numbers of cells
numCellsSemiInfinite_pos = ceil(6);
numCellsSemiInfinite_neg = ceil(6);

numCellsInfinite_pos = ceil(6*cutvecpos(1));
numCellsInfinite_neg = ceil(6*cutvecneg(1));

numFloquetPoints_pos = 40;
numFloquetPoints_neg = 40;
numFloquetPoints_pos = max(numFloquetPoints_pos, ceil(numFloquetPoints_pos/cutvecpos(1)));
numFloquetPoints_neg = max(numFloquetPoints_neg, ceil(numFloquetPoints_neg/cutvecneg(1)));

FloquetPoints_pos = linspace(-pi, pi, numFloquetPoints_pos).';
FloquetPoints_neg = linspace(-pi, pi, numFloquetPoints_neg).';

%% Meshes
struct_mesh = 0;

numNodesX = 40;
numNodesYpos = 40; numNodesYpos = max(numNodesYpos, ceil(numNodesYpos/cutvecpos(1)));
numNodesYneg = 40; numNodesYneg = max(numNodesYneg, ceil(numNodesYneg/cutvecneg(1)));
numNodesZ = 40;
numNodes_int = numNodesZ/4;
opts.ratio_quad = 4;

meshXYpos = meshes.MeshRectangle(struct_mesh, [0  1], [0 1/cutvecpos(1)], numNodesX, numNodesYpos);
meshXYneg = meshes.MeshRectangle(struct_mesh, [0 -1], [0 1/cutvecneg(1)], numNodesX, numNodesYneg);

meshZXpos = meshes.MeshRectangle(struct_mesh, [0 1], [0  1], numNodesZ, numNodesX);
meshZXneg = meshes.MeshRectangle(struct_mesh, [0 1], [0 -1], numNodesZ, numNodesX);

meshYZpos = meshes.MeshRectangle(struct_mesh, [0 1], [0  1], numNodesYpos, numNodesZ);
meshYZneg = meshes.MeshRectangle(struct_mesh, [0 1], [0  1], numNodesYneg, numNodesZ);

meshLineZ = meshes.MeshSegment('uniform', 0, 1, numNodesZ);
meshInter = meshes.MeshSegment('uniform', -Lint, Lint, floor(2*Lint*numNodes_int));

%% Basis functions and boundary conditions
basis_functions_X = 'Fourier';
basis_functions_Y = 'Fourier';

% Positive side
if strcmpi(basis_functions_X, 'Lagrange')
  BCstruct.pos.spBX = FEPack.spaces.PeriodicLagrangeBasis(meshYZpos.domain('volumic')); % X = cst
else
  BCstruct.pos.spBX = FEPack.spaces.FourierBasis(meshYZpos.domain('volumic'), [floor(numNodesYpos/4), floor(numNodesZ/4)]);
end

if strcmpi(basis_functions_Y, 'Lagrange')
  % Particular Lagrange for basis functons on Y = cst
  ecs_fun = @(u, dom) ((((u|dom.mesh.domain('xmin')) - (u|dom.mesh.domain('xmax'))) == 0.0) &...      % Periodic on Z
                        ((u|dom.mesh.domain('ymin')) == 0.0) & ((u|dom.mesh.domain('ymax')) == 0.0)); % Dirichlet on X
  
  BCstruct.pos.spBY = FEPack.spaces.EssentialBCLagrangeBasis(meshZXpos.domain('volumic'), ecs_fun); % Y = cst
else
  BCstruct.pos.spBY = FEPack.spaces.FourierBasis(meshZXpos.domain('volumic'), [floor(numNodesZ/4), floor(numNodesX/4)]); % Y = cst
end

% Other parameters
BCstruct.pos.BCdu = 0.0;
BCstruct.pos.BCu = 1.0;
% BCstruct.pos.representation = 'projection';

% Negative side
BCstruct.neg = BCstruct.pos;
if strcmpi(basis_functions_X, 'Lagrange')
  BCstruct.neg.spBX = FEPack.spaces.PeriodicLagrangeBasis(meshYZneg.domain('volumic')); % X = cst
else
  BCstruct.neg.spBX = FEPack.spaces.FourierBasis(meshYZneg.domain('volumic'), [floor(numNodesYneg/4), floor(numNodesZ/4)]);
end

if strcmpi(basis_functions_Y, 'Lagrange')
  % Particular Lagrange for basis functons on Y = cst
  ecs_fun = @(u, dom) ((((u|dom.mesh.domain('xmin')) - (u|dom.mesh.domain('xmax'))) == 0.0) &...      % Periodic on Z
                        ((u|dom.mesh.domain('ymin')) == 0.0) & ((u|dom.mesh.domain('ymax')) == 0.0)); % Dirichlet on X

  BCstruct.neg.spBY = FEPack.spaces.EssentialBCLagrangeBasis(meshZXneg.domain('volumic'), ecs_fun); % Y = cst
else
  BCstruct.neg.spBY = FEPack.spaces.FourierBasis(meshZXneg.domain('volumic'), [floor(numNodesZ/4), floor(numNodesX/4)]); % Y = cst
end



% close all;
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');

% mu2Dglob =  @(x) (x(:, 1) >= 0) .*  mu2Dpos(x) + (x(:, 1) <  0) .* mu2Dneg(x);
% rho2Dglob = @(x) (x(:, 1) >= 0) .* rho2Dpos(x) + (x(:, 1) <  0) .* rho2Dneg(x);

% % opts.numCellsInfinite = 3;

% mesh = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0 opts.period], numNodesX, numNodesX);

% H = figure('Position', get(0, 'Screensize'), 'visible', 'off');
% for idS = 1:(2*numCellsSemiInfinite_pos)
%   for idI = 1:(2*numCellsInfinite_pos)
%     X = mesh.points(:, 1) + (idS - numCellsSemiInfinite_pos - 1);
%     Y = mesh.points(:, 2) + (idI - numCellsInfinite_pos - 1);

%     trisurf(mesh.triangles, X, Y, rho2Dglob([X, Y]));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%   end
% end
% view(2)
% axis off
% colorbar off;
% axis([-6 6 -6 6]);
% caxis([0.5, 1.5]);
% print('last_simus/medium', '-dpng');

% error();













%% Compute DtN operators
% //////////////////////
opts.suffix = 'pos';
opts.cutvec = cutvecpos;
opts.cutmat = cutmat_pos;
opts.numFloquetPoints = numFloquetPoints_pos;
opts.FloquetPoints = FloquetPoints_pos;
quasi2DHalfGuideDirichlet_par(+1, meshXYpos, meshLineZ, mu3Dpos, rho3Dpos, BCstruct.pos, opts);
system(['rm ', nomdossier, 'FEmat_*']);
%%
opts.suffix = 'neg';
opts.cutvec = cutvecneg;
opts.cutmat = cutmat_neg;
opts.numFloquetPoints = numFloquetPoints_neg;
opts.FloquetPoints = FloquetPoints_neg;
quasi2DHalfGuideDirichlet_par(-1, meshXYneg, meshLineZ, mu3Dneg, rho3Dneg, BCstruct.neg, opts);
system(['rm ', nomdossier, 'FEmat_*']);

%% Basis functions associated to the interface
% ////////////////////////////////////////////
type_basis_int = 'sine'; % 'sine', 'Lagrange'
type_basis_trs = 'Fourier'; % 'Fourier'

% Interface direction
if strcmpi(type_basis_int, 'sine')

  phisInt = @(x, n) sqrt(1 / Lint) * sin(pi * (x(:, 1) + Lint) * n / (2 * Lint)) .* (abs(x(:, 1)) <= Lint);
  numBasisInt = floor(meshInter.numPoints / 2);

elseif strcmpi(type_basis_int, 'Lagrange')

  % error('Il n''y a pas de bug. Je veux juste que tu vérifies le nombre de noeuds (pas trop grand).');
  phisInt = @(x, n) ((x(:, 1) - meshInter.points(n, 1).') ./ (meshInter.points(n+1, 1).' - meshInter.points(n, 1).'))...
                                  .* (x(:, 1) >= meshInter.points(n, 1).' & x(:, 1) < meshInter.points(n+1, 1).')...
                               + ((meshInter.points(n+2, 1).' - x(:, 1)) ./ (meshInter.points(n+2, 1).' - meshInter.points(n+1, 1).'))...
                                  .* (x(:, 1) >= meshInter.points(n+1, 1).' & x(:, 1) < meshInter.points(n+2, 1).');
  numBasisInt = meshInter.numPoints - 2;

else

  error(['Type de fonction de base', type_basis_int, 'non reconnu.']);

end
spBint = FEPack.spaces.SpectralBasis(meshInter.domain('volumic'), phisInt, numBasisInt);
spBint.computeBasisMatrices;

% Transverse direction
if strcmpi(type_basis_trs, 'Fourier')
  spBtrs = FEPack.spaces.FourierBasis(meshLineZ.domain('volumic'), numNodesZ/4);
else
  error(['Type de fonction de base', type_basis_trs, 'non reconnu.']);
end

numBasis = spBint.numBasis * spBtrs.numBasis;
[idInt, idTrs] = ind2sub([spBint.numBasis, spBtrs.numBasis], 1:numBasis);

phisYZpos = @(YZ, idI) spBtrs.phis(YZ(:, 2) - YZ(:, 1) * cutslopePos, idTrs(idI)) .* phisInt(YZ(:, 1) / cutvecpos(1), idInt(idI));
phisYZneg = @(YZ, idI) spBtrs.phis(YZ(:, 2) - YZ(:, 1) * cutslopeNeg, idTrs(idI)) .* phisInt(YZ(:, 1) / cutvecneg(1), idInt(idI));

%% Positive half-space DtN
% ////////////////////////
fprintf('<strong>DtN demi-plan positif.</strong>\n');
TFBphiPosVec = zeros(meshYZpos.numPoints, numBasis, numFloquetPoints_pos);
for idI = 1:numBasis
  TFBphiPosVec(:, idI, :) = FEPack.tools.BlochTransform(meshYZpos.points,...
                                           FloquetPoints_pos, @(x) phisYZpos(x, idI),...
                                           1, opts.period, 1000);
end

TFBphiPos = cell(numFloquetPoints_pos, 1);
lambda_pos = zeros(numBasis);
wpos = (2*pi / opts.period) / (numFloquetPoints_pos - 1);
FE_to_spectral = BCstruct.pos.spBX.FE_to_spectral;
massmat = BCstruct.pos.spBX.massmat;

parfor idFB = 1:numFloquetPoints_pos
  fprintf('%d sur %d\n', idFB, numFloquetPoints_pos);
  solguidePos = load([nomdossier, 'half_guide_sol_pos_Floquet_', num2str(idFB)]);

  TFBphiPos{idFB} = FE_to_spectral * TFBphiPosVec(:, :, idFB);
  TFBlambdaPos = massmat * solguidePos.Lambda;

  lambda_pos = lambda_pos + wpos * TFBphiPos{idFB}' * TFBlambdaPos * TFBphiPos{idFB};
end

lambda_pos = lambda_pos / cutvecpos(1);

%% Negative half-space DtN
% ////////////////////////
fprintf('<strong>DtN demi-plan negatif.</strong>\n');
TFBphiNegVec = zeros(meshYZneg.numPoints, numBasis, numFloquetPoints_neg);
for idI = 1:numBasis
  TFBphiNegVec(:, idI, :) = FEPack.tools.BlochTransform(meshYZneg.points,...
                                           FloquetPoints_neg, @(x) phisYZneg(x, idI),...
                                           1, opts.period, 1000);
end

TFBphiNeg = cell(numFloquetPoints_neg, 1);
lambda_neg = zeros(numBasis);
wneg = (2*pi / opts.period) / (numFloquetPoints_neg - 1);
FE_to_spectral = BCstruct.neg.spBX.FE_to_spectral;
massmat = BCstruct.neg.spBX.massmat;

parfor idFB = 1:numFloquetPoints_neg
  fprintf('%d sur %d\n', idFB, numFloquetPoints_neg);
  solguideNeg = load([nomdossier, 'half_guide_sol_neg_Floquet_', num2str(idFB)]);

  TFBphiNeg{idFB} = FE_to_spectral * TFBphiNegVec(:, :, idFB);
  TFBlambdaNeg = massmat * solguideNeg.Lambda;

  lambda_neg = lambda_neg + wneg * TFBphiNeg{idFB}' * TFBlambdaNeg * TFBphiNeg{idFB};
end

lambda_neg = lambda_neg / cutvecneg(1);

%% Solve the integral equation on the interface
% /////////////////////////////////////////////
fprintf('<strong>Equation d''interface.</strong>\n');
numNodesEF = meshInter.numPoints * meshLineZ.numPoints;
[idIntEF, idTrsEF] = ind2sub([meshInter.numPoints, meshLineZ.numPoints], 1:numNodesEF);

% massmat = spBint.massmat(idInt, idInt) .* spBtrs.massmat(idTrs, idTrs);
% lambda_pos = massmat; lambda_neg = massmat;

projmat = spBint.projmat(idInt, idIntEF) .* spBtrs.projmat(idTrs, idTrsEF);
Gvec = Gint(meshInter.points(idIntEF, :)) .* Gtrs(meshLineZ.points(idTrsEF, :));
mat_G_v = projmat * Gvec;

% -------------------------------------------------------------------- %
% The minus sign comes from the definition of the Lambda when they are %
% the DtN operators (see PeriodicHalfGuideBVP.m)                       %
% -------------------------------------------------------------------- %
trace_solution = -(lambda_pos + lambda_neg) \ mat_G_v;                 %
% -------------------------------------------------------------------- %

%% FB transform of trace
% //////////////////////
parfor idFB = 1:numFloquetPoints_pos
  fprintf('%d sur %d\n', idFB, numFloquetPoints_pos);
  
  soltrace = [];
  soltrace.vec = TFBphiPos{idFB} * trace_solution;
  
  % Save the trace of the solution
  FEPack.tools.parsave([nomdossier, 'sol_trace_Floquet_pos_', num2str(idFB)], soltrace, true);
end

parfor idFB = 1:numFloquetPoints_neg
  fprintf('%d sur %d\n', idFB, numFloquetPoints_neg);
  
  soltrace = [];
  soltrace.vec = TFBphiNeg{idFB} * trace_solution;
  
  % Save the trace of the solution
  FEPack.tools.parsave([nomdossier, 'sol_trace_Floquet_neg_', num2str(idFB)], soltrace, true);
end

%% Construct solution
opts.suffix = 'pos';
opts.cutvec = cutvecpos;
opts.cutmat = cutmat_pos;
opts.numFloquetPoints = numFloquetPoints_pos;
opts.FloquetPoints = FloquetPoints_pos;
opts.numCellsInfinite = numCellsInfinite_pos;
opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
construct_solution_par(meshXYpos, meshLineZ, BCstruct.pos, opts, 'sol_trace_Floquet_pos_');

%
opts.suffix = 'neg';
opts.cutvec = cutvecneg;
opts.cutmat = cutmat_neg;
opts.numFloquetPoints = numFloquetPoints_neg;
opts.FloquetPoints = FloquetPoints_neg;
opts.numCellsInfinite = numCellsInfinite_neg;
opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
construct_solution_par(meshXYneg, meshLineZ, BCstruct.neg, opts, 'sol_trace_Floquet_neg_');

%% Plot solution
mafig = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idY = 1:2*numCellsInfinite_pos
  sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);

  for idX = 1:numCellsSemiInfinite_pos
    X = meshXYpos.points(:, 1) + (idX - 1);
    Y = meshXYpos.points(:, 2) + (idY - numCellsInfinite_pos - 1) / cutvecpos(1);

    trisurf(meshXYpos.triangles, X, Y, real(sol.Usol(:, idX)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end
%
%
for idY = 1:2*numCellsInfinite_neg
  sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);

  for idX = 1:numCellsSemiInfinite_neg
    X = meshXYneg.points(:, 1) - (idX - 1);
    Y = meshXYneg.points(:, 2) + (idY - numCellsInfinite_neg - 1) / cutvecneg(1);

    trisurf(meshXYneg.triangles, X, Y, real(sol.Usol(:, idX)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

% xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
% ylim([-opts.numCellsInfinite/cutvecpos(1), opts.numCellsInfinite/cutvecneg(1)]);

axis([-6 6 -6 6]);

% caxis([-0.25, 0.25]);

% system(['rm ', nomdossier, 'half_guide_*']);
% system(['rm ', nomdossier, 'local_cell_sol_*']);
% system(['rm ', nomdossier, 'sol_trace_Floquet_*']);

savefig(mafig, ['last_simus/solution_', int2str(idtest)], 'compact');
% print('last_simus/solution', '-dpng')
end