%%
problem_setting  = entrees{1};
period           = entrees{2};
numNodes         = entrees{3};
test_omega       = entrees{4};
test_case        = entrees{5};
test_per         = entrees{6};
test_G           = entrees{7};
compute_3D_sol   = entrees{8};
compare_w_refsol = entrees{9};
delete_sol       = entrees{10};
suffix           = entrees{11};
theta            = entrees{12};
  
import FEPack.*
profile OFF

if (length(test_per) == 1)
  test_per = [1, test_per];
end

% opts.nomdossier = ['last_simus/solution_', problem_setting, '_', test_case, '_', test_omega, '_per', num2str(test_per), '_', test_G, suffix, '/']; % /!\ Pas oublier le / à la fin
opts.nomdossier = ['/UMA/tmp/amenoagbadji/solution_', problem_setting, '_', test_case, '_', test_omega, '_per', num2str(test_per(1)), num2str(test_per(2)), '_', test_G, '/']; % /!\ Pas oublier le / à la fin
nomdossier = opts.nomdossier;

system(['mkdir ', nomdossier]);

% if ~(compare_w_refsol)
%   warning(['Fichiers dans ', nomdossier, ' supprimés.']);
%   system(['rm ', nomdossier, '*']);
% end

%% Problem-related variables
opts.period = period;
opts.verbose = 0;

if strcmpi(test_omega, 'omega1')
  opts.omega = 8 + 0.25i;
elseif strcmpi(test_omega, 'omega2')
  opts.omega = 20 + 0.25i;
elseif strcmpi(test_omega, 'omega3')
  opts.omega = 20 + 0.05i;
elseif strcmpi(test_omega, 'omega4')
  opts.omega = 1 + 0.25i;
end

mu2Dpos =  @(x) ones(size(x, 1), 1);
mu2Dneg  = @(x) ones(size(x, 1), 1);

if strcmpi(problem_setting, 'A')

  % 2D coefficients
  period_posFun = 1;

  if strcmpi(test_case, 'rational')
    period_negFun = 1;
  else
    period_negFun = sqrt(2);
  end

  % echelle = 0.1;
  if strcmpi(test_case, 'constant')
    rho2Dpos = @(x) ones(size(x, 1), 1);
    rho2Dneg = @(x) 2 * ones(size(x, 1), 1);
  else
    rho2DposCell = @(x) 0.5 + tools.FUN_per_cutoff_circle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.4, 0.4]);
    rho2DnegCell = @(x) 0.5 + tools.FUN_per_cutoff_cuboid(x, [1; 0], [0; period_negFun], [0.5, 0.5], [-0.25, 0.25], [-0.25, 0.25], 0.5, 1);
    
    rho2Dpos = @(x) rho2DposCell(x / opts.period);
    rho2Dneg = @(x) rho2DnegCell(x / opts.period);
  end

  % Cut vector
  test_per
  period_pos = test_per(1) * period_posFun;
  period_neg = test_per(2) * period_negFun;

  opts.cutvec = [1/period_pos; 1/period_neg];

  % 3D coefficients
  mu3Dpos  = @(x) mu2Dpos( [x(:, 1), period_pos * x(:, 2)]);
  mu3Dneg  = @(x) mu2Dneg( [x(:, 1), period_neg * x(:, 3)]);
  rho3Dpos = @(x) rho2Dpos([x(:, 1), period_pos * x(:, 2)]);
  rho3Dneg = @(x) rho2Dneg([x(:, 1), period_neg * x(:, 3)]);
  
else

  % 2D coefficients
  if strcmpi(test_case, 'rational')
    vecperFun = [1, 1];
  else
    vecperFun = [cos(3*pi/5), sin(3*pi/5)];
  end
  
  v1 = [ cos(theta); sin(theta)];
  v2 = [-sin(theta); cos(theta)];

  if (abs(v1(1)) > abs(v2(1)))
    vI = [ abs(v1(1)); abs(v1(2))];
    vO = [-abs(v2(1)); abs(v2(2))];
  else
    vI = [ abs(v2(1)); abs(v2(2))];
    vO = [-abs(v1(1)); abs(v1(2))];
  end

  % % theta entre 0 et 2*pi (exclus);
  % if (theta < pi/4) || (theta >= 3*pi/4 && theta < 5*pi/4) || (theta >= 7*pi/4)
  %   vI = [ abs(v1(1)); abs(v1(2))];
  %   vO = [-abs(v2(1)); abs(v2(2))];
  % else
  %   vI = [ abs(v2(1)); abs(v2(2))];
  %   vO = [-abs(v1(1)); abs(v1(2))];
  % end

  % sgnvIx = sign(vI(1));
  % sgnvPz = sign(vO(2) - (vO(1) * vI(2) / vI(1)));
  % invT = [sgnvIx * vI, [0; 1] * sgnvPz];
  % vI = v1;
  % vO = v2;

  invT = [vI, [0; 1]];
  T = inv(invT);
  vecperFun = abs(T * vO);

  invRotmat = inv([v1, v2]);

  if strcmpi(test_case, 'constant')
    rho2Dpos = @(x) ones(size(x, 1), 1);
    rho2Dneg = @(x) 2 * ones(size(x, 1), 1);
  else
    % rho2DposInit = @(x) 0.5 + tools.FUN_per_cutoff_circle((rotmat * x')'/opts.period, [1; 0], vecperFun, [0.5, 0.5], [-0.4, 0.4]);
    % rho2DposInit = @(x) 0.5 + tools.FUN_per_cutoff_circle((invRotmat * x')'/opts.period, [1; 0], [0; 1], [0.4, 0.6], [-0.4, 0.4]);

    rho2DposInit = @(x) 0.5 + tools.FUN_per_cutoff_kite((invRotmat * x')'/opts.period, [1; 0], [0; 1], [0.5, 0.5], 1.5*[1, 1.3, 1.5]/6);
    rho2Dneg = @(x) ones(size(x, 1), 1);

    rho2Dpos = @(x) rho2DposInit((invT * x')');
    mu2Dpos =  @(x) kron(T * T', ones(size(x, 1), 1));
    mu2Dneg =  @(x) kron(T * T', ones(size(x, 1), 1));
  end

  BBx = 6; BBy = 6;
  nint = ceil(BBx / abs(vI(1)));
  mint = ceil((BBy + nint * abs(vI(2))) / abs(vecperFun(2)));
  
  % Cut vector
  % vecper = test_per(2) * [-vecperFun(1), vecperFun(2)];
  vecper = [vecperFun(1), vecperFun(2)];
  opts.cutvec = [1/vecper(2); -vecper(1)/vecper(2)];
  
  % 3D coefficients
  fun3D = @(fun2D, x) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
  mu3Dpos  = @(x) fun3D(mu2Dpos,  x);
  mu3Dneg  = @(x) fun3D(mu2Dneg,  x);
  rho3Dpos = @(x) fun3D(rho2Dpos, x);
  rho3Dneg = @(x) fun3D(rho2Dneg, x);

end

% Cut matrix and cut slope
opts.cutmat = [[1; 0; 0], [0; opts.cutvec]];

% Jump data
G = @(x) FEPack.tools.cutoff(x, -0.5, 0.5);

if strcmpi(problem_setting, 'B')
  G = @(x) G(x) * vI(1);
end

if strcmpi(test_G, 'G1')
  G3D = @(x) G(x(:, 1)/opts.cutvec(1));
else
  G3D = @(x) G(x(:, 1)/opts.cutvec(1)) .* exp(2i*pi*(x(:, 2) - x(:, 1)*opts.cutvec(2)/opts.cutvec(1)));
end

% Numbers of cells
numCellsSemiInfinite_pos = ceil(nint/opts.period);
numCellsSemiInfinite_neg = ceil(nint/opts.period);
opts.numCellsInfinite = ceil(mint/opts.period);
opts.numFloquetPoints = ceil(60/opts.period);
opts.numFloquetPoints = 2*ceil(opts.numFloquetPoints/2);
opts.numFloquetPoints = max(opts.numFloquetPoints, ceil(opts.numFloquetPoints/opts.cutvec(1)));

if (opts.numFloquetPoints == 1)
  opts.FloquetPoints = 0;
else
  opts.FloquetPoints = linspace(-pi/period, pi/period, opts.numFloquetPoints+1);
  opts.FloquetPoints(end) = [];
end

%% Meshes
struct_mesh = 1;
numNodesX = numNodes;
numNodesY = numNodes;
numNodesZ = numNodes;
opts.ratio_quad = 8;

meshXYpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period/opts.cutvec(1)], numNodesX, max(numNodesY, ceil(numNodesY/abs(opts.cutvec(1)))));
meshXYneg = FEPack.meshes.MeshRectangle(struct_mesh, [0 -opts.period], [0  opts.period/opts.cutvec(1)], numNodesX, max(numNodesY, ceil(numNodesY/abs(opts.cutvec(1)))));
meshYZ    = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesY, numNodesZ);
meshLineZ = FEPack.meshes.MeshSegment('uniform', 0, opts.period, numNodesZ);

% %% Plot the coefficients
% mafig = figure;
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');

% for idY = 1:2*opts.numCellsInfinite
%   for idX = 1:numCellsSemiInfinite_pos
%     X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
%     Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

%     trisurf(meshXYpos.triangles, X, Y, rho2Dpos([X, Y]));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%   end
% end

% for idY = 1:2*opts.numCellsInfinite
%   for idX = 1:numCellsSemiInfinite_neg
%     X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
%     Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

%     trisurf(meshXYneg.triangles, X, Y, rho2Dneg([X, Y]));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%   end
% end

% savefig(mafig, [nomdossier, 'fig_coeff_2D'], 'compact');
% close(mafig);

%% Basis functions and boundary conditions
% Positive side
basis_functions_X = 'Fourier';
basis_functions_S = 'Fourier';

if strcmpi(basis_functions_X, 'Lagrange')
  BCstruct.pos.spBX = FEPack.spaces.PeriodicLagrangeBasis(meshYZ.domain('volumic')); % X = cst
else
  nBY = max(1, floor(numNodesY/4));
  nBZ = max(1, floor(numNodesZ/4));
  BCstruct.pos.spBX = FEPack.spaces.FourierBasis(meshYZ.domain('volumic'), [nBY, nBZ]);
end

if strcmpi(basis_functions_S, 'Lagrange')
  BCstruct.pos.spBS = FEPack.spaces.PeriodicLagrangeBasis(meshLineZ.domain('volumic')); % S/Z line
else
  nBS = max(1, floor(numNodesZ/4));
  BCstruct.pos.spBS = FEPack.spaces.FourierBasis(meshLineZ.domain('volumic'), nBS);
end

% Other parameters
BCstruct.pos.BCdu = 0.0;
BCstruct.pos.BCu = 1.0;

% Negative side
BCstruct.neg = BCstruct.pos;

%% Compute DtN operators
opts.suffix = 'pos';
quasi2DHalfGuideDirichlet_par(+1, meshXYpos, meshLineZ, mu3Dpos, rho3Dpos, BCstruct.pos, opts);
system(['rm ', nomdossier, 'FEmat_*']);
%
opts.suffix = 'neg';
quasi2DHalfGuideDirichlet_par(-1, meshXYneg, meshLineZ, mu3Dneg, rho3Dneg, BCstruct.neg, opts);
system(['rm ', nomdossier, 'FEmat_*']);

%% Interface equation
fprintf('<strong>Equation d''interface.</strong>\n');
FloquetPoints = opts.FloquetPoints;
numFloquetPoints = opts.numFloquetPoints;
pointsYZ = meshYZ.points;
spBX_FE_to_spectral = BCstruct.pos.spBX.FE_to_spectral;
period = opts.period;

parfor idFB = 1:numFloquetPoints
  fprintf('%d sur %d\n', idFB, numFloquetPoints);
  FloquetVar = FloquetPoints(idFB);

  % DtN operators
  solguidePos = load([nomdossier, 'half_guide_sol_pos_Floquet_', num2str(idFB)]);
  solguideNeg = load([nomdossier, 'half_guide_sol_neg_Floquet_', num2str(idFB)]);
  
  Lambda_pos = solguidePos.Lambda;
  Lambda_neg = solguideNeg.Lambda;

  % Rhs
  jumpData_FB = @(x) FEPack.tools.BlochTransform(x, FloquetVar, G3D, 1, period, 1000);
  GG = spBX_FE_to_spectral * jumpData_FB(pointsYZ);

  % The minus sign comes from the definition of the Lambda
  soltrace = [];
  soltrace.vec = -(Lambda_pos + Lambda_neg) \ GG;

  % Save the trace of the solution
  FEPack.tools.parsave([nomdossier, 'sol_trace_Floquet_', num2str(idFB)], soltrace, true);
end

%% Construct solution
opts.suffix = 'pos';
opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
construct_solution_par(+1, meshXYpos, meshLineZ, BCstruct.pos, opts);

opts.suffix = 'neg';
opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
construct_solution_par(-1, meshXYneg, meshLineZ, BCstruct.neg, opts);

% %% Plot solution
if strcmpi(test_case, 'homogeneisation')
  close all;
  numCellsInfinite = opts.numCellsInfinite;

  % Compute limit solution
  % //////////////////////
  xibound = 500;
  xi = linspace(-xibound, xibound, floor(10 * (2 * xibound))).';
  Delta_xi = xi(2) - xi(1);

  % Fourier transform of interface data
  ypts = linspace(-0.5, 0.5, ceil(10*xibound)).';
  Delta_y = ypts(2) - ypts(1);
  FT_G = Delta_y * exp(-1i * xi * ypts.') * G(ypts);
  
  % //////////////////// %
  % Configuration A only %
  % //////////////////// %
  rho_eff_pos = integral2(@(X, Y) reshape(rho2DposCell([X(:), Y(:)]), size(X)), 0, 1, 0, period_posFun) / period_posFun;
  rho_eff_neg = integral2(@(X, Y) reshape(rho2DnegCell([X(:), Y(:)]), size(X)), 0, 1, 0, period_negFun) / period_negFun;
  %
  coe_pos = sqrt(xi.^2 - rho_eff_pos * opts.omega^2);
  coe_neg = sqrt(xi.^2 - rho_eff_neg * opts.omega^2);
  %
  solrefPos = @(X, Y) (Delta_xi/(2*pi)) * exp(-X * coe_pos.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));
  solrefNeg = @(X, Y) (Delta_xi/(2*pi)) * exp(+X * coe_neg.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));
  %
  %
  %
  mafig = figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  for idY = 1:2*opts.numCellsInfinite
    sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);

    for idX = 1:numCellsSemiInfinite_pos
      X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
      Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

      trisurf(meshXYpos.triangles, X, Y, real(sol.Usol(:, idX)));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end
  %
  %
  for idY = 1:2*opts.numCellsInfinite
    sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);

    for idX = 1:numCellsSemiInfinite_neg
      X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
      Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

      trisurf(meshXYneg.triangles, X, Y, real(sol.Usol(:, idX)));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end

  xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
  ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
  savefig(mafig, [nomdossier, 'fig_solution_2D'], 'compact');
  close(mafig);


  mafig = figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  for idY = 1:2*opts.numCellsInfinite
    sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);

    for idX = 1:numCellsSemiInfinite_pos
      X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
      Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

      solErreur = sol.Usol(:, idX) - solrefPos(X, Y);
      trisurf(meshXYpos.triangles, X, Y, real(solErreur));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end
  %
  %
  for idY = 1:2*opts.numCellsInfinite
    sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);

    for idX = 1:numCellsSemiInfinite_neg
      X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
      Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);

      solErreur = sol.Usol(:, idX) - solrefNeg(X, Y);
      trisurf(meshXYneg.triangles, X, Y, real(solErreur));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end

  xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
  ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
  savefig(mafig, [nomdossier, 'fig_erreur_2D'], 'compact');
  close(mafig);

elseif (strcmpi(test_case, 'constant') & compare_w_refsol)
  % Compute reference solution
  % //////////////////////
  xibound = 500;
  xi = linspace(-xibound, xibound, floor(10 * (2 * xibound))).';
  Delta_xi = xi(2) - xi(1);

  % Fourier transform of interface data
  ypts = linspace(-0.5, 0.5, ceil(10*xibound)).';
  Delta_y = ypts(2) - ypts(1);
  FT_G = Delta_y * exp(-1i * xi * ypts.') * G(ypts);

  % //////////////////// %
  % Configuration A only %
  % //////////////////// %
  rho_const_pos = rho2Dpos([0, 0]);
  rho_const_neg = rho2Dneg([0, 0]);
  %
  coe_pos = sqrt(xi.^2 - rho_const_pos * opts.omega^2);
  coe_neg = sqrt(xi.^2 - rho_const_neg * opts.omega^2);
  %
  solrefPos = @(X, Y) (Delta_xi/(2*pi)) * exp(-X * coe_pos.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));
  solrefNeg = @(X, Y) (Delta_xi/(2*pi)) * exp(+X * coe_neg.' + 1i * Y * xi.') * (FT_G ./ (coe_pos + coe_neg));

  u = FEPack.pdes.PDEObject; v = dual(u);

  mat_gradu_gradv_pos = FEPack.pdes.Form.intg(meshXYpos.domain('volumic'), grad2(u) * grad2(v));
  mat_u_v_pos         = FEPack.pdes.Form.intg(meshXYpos.domain('volumic'), u * v);
  mat_gradu_gradv_neg = FEPack.pdes.Form.intg(meshXYneg.domain('volumic'), grad2(u) * grad2(v));
  mat_u_v_neg         = FEPack.pdes.Form.intg(meshXYneg.domain('volumic'), u * v);
  
  mafig = figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  ecartL2pos = 0; ecartH1pos = 0;
  normUL2pos = 0; normUH1pos = 0;
  for idY = 1:2*opts.numCellsInfinite
    sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);
    for idX = 1:numCellsSemiInfinite_pos
      X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
      Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);
      ecart = sol.Usol(:, idX) - solrefPos(X, Y);

      ecartL2pos = ecartL2pos + (ecart' * mat_u_v_pos * ecart);
      normUL2pos = normUL2pos + (solrefPos(X, Y)' * mat_u_v_pos * solrefPos(X, Y));

      ecartH1pos = ecartH1pos + (ecart' * mat_gradu_gradv_pos * ecart);
      normUH1pos = normUH1pos + (solrefPos(X, Y)' * mat_gradu_gradv_pos * solrefPos(X, Y));

      trisurf(meshXYpos.triangles, X, Y, real(ecart));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end
  fid = fopen([nomdossier, 'erreurs_pos.txt'], 'a+');
  fprintf(fid, '%d\t%0.5d\t%0.5d\t%0.5d\t%0.5d\n', numNodesX, ecartL2pos, normUL2pos, ecartH1pos, normUH1pos);
  fclose(fid);

  ecartL2neg = 0; ecartH1neg = 0;
  normUL2neg = 0; normUH1neg = 0;
  for idY = 1:2*opts.numCellsInfinite
    sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);
    for idX = 1:numCellsSemiInfinite_neg
      X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
      Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);
      ecart = sol.Usol(:, idX) - solrefNeg(X, Y);

      ecartL2neg = ecartL2neg + (ecart' * mat_u_v_neg * ecart);
      normUL2neg = normUL2neg + (solrefNeg(X, Y)' * mat_u_v_neg * solrefNeg(X, Y));

      ecartH1neg = ecartH1neg + (ecart' * mat_gradu_gradv_neg * ecart);
      normUH1neg = normUH1neg + (solrefNeg(X, Y)' * mat_gradu_gradv_neg * solrefNeg(X, Y));

      trisurf(meshXYneg.triangles, X, Y, real(ecart));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end
  fid = fopen([nomdossier, 'erreurs_neg.txt'], 'a+');
  fprintf(fid, '%d\t%0.5d\t%0.5d\t%0.5d\t%0.5d\n', numNodesX, ecartL2neg, normUL2neg, ecartH1neg, normUH1neg);
  fclose(fid);

  xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
  ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
  % caxis([-0.25, 0.25]);
  savefig(mafig, [nomdossier, 'ecart_solution_2D'], 'compact');

else
  mafig = figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  theNs = (1:nint);
  theMs = (-mint:mint-1);
  X = ones(2*mint, 1) * (theNs - 1);
  Y = theMs' * ones(1, nint) / opts.cutvec(1);
  Pnew = (invT * [X(:), Y(:)]')';
  Xnew = reshape(Pnew(:, 1), 2*mint, nint);
  Ynew = reshape(Pnew(:, 2), 2*mint, nint);

  Iin_ind = find(...% (Xnew + abs(vI(1)) >= -BBx) &...
                 (Xnew - abs(vI(1)) <= BBx) &...
                 (Ynew + abs(vI(2)) +  abs(vecperFun(2)) >= -BBy) &...
                 (Ynew - abs(vI(2)) -  abs(vecperFun(2)) <=  BBy));
  [IdY, IdX] = ind2sub([2*mint, nint], Iin_ind);

  for idI = 1:length(Iin_ind)
    sol = load([nomdossier, 'solution_pos_Y_', num2str(theMs(IdY(idI)) + mint + 1)]);

    X = meshXYpos.points(:, 1) + (theNs(IdX(idI)) - 1) * opts.period;
    Y = meshXYpos.points(:, 2) + (theMs(IdY(idI))) * opts.period / opts.cutvec(1);

    Pnew = (invT * [X, Y]')';
    Xnew = Pnew(:, 1);
    Ynew = Pnew(:, 2);

    trisurf(meshXYpos.triangles, Xnew, Ynew, real(sol.Usol(:, theNs(IdX(idI)))));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end

  theNs = (1:nint);
  theMs = (-mint:mint-1);
  X = ones(2*mint, 1) * -1 * (theNs - 1);
  Y = theMs' * ones(1, nint) / opts.cutvec(1);
  Pnew = (invT * [X(:), Y(:)]')';
  Xnew = reshape(Pnew(:, 1), 2*mint, nint);
  Ynew = reshape(Pnew(:, 2), 2*mint, nint);

  Iin_ind = find((Xnew + abs(vI(1)) >= -BBx) &...
                 ...%  (Xnew - abs(vI(1)) <=  BBx) &...
                 (Ynew + abs(vI(2)) + abs(1/opts.cutvec(1)) >= -BBy) &...
                 (Ynew - abs(vI(2)) - abs(1/opts.cutvec(1)) <=  BBy));
  [IdY, IdX] = ind2sub([2*mint, nint], Iin_ind);

  for idI = 1:length(Iin_ind)
    sol = load([nomdossier, 'solution_neg_Y_', num2str(theMs(IdY(idI)) + mint + 1)]);

    X = meshXYneg.points(:, 1) - (theNs(IdX(idI)) - 1) * opts.period;
    Y = meshXYneg.points(:, 2) + (theMs(IdY(idI))) * opts.period / opts.cutvec(1);

    Pnew = (invT * [X, Y]')';
    Xnew = Pnew(:, 1);
    Ynew = Pnew(:, 2);

    trisurf(meshXYneg.triangles, Xnew, Ynew, real(sol.Usol(:, theNs(IdX(idI)))));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end

  % for idY = 1:2*opts.numCellsInfinite
  %   sol = load([nomdossier, 'solution_pos_Y_', num2str(idY)]);

  %   for idX = 1:numCellsSemiInfinite_pos
  %     X = meshXYpos.points(:, 1) + (idX - 1) * opts.period;
  %     Y = meshXYpos.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);
  %     Pnew = (invT * [X, Y]')';
  %     Xnew = Pnew(:, 1);
  %     Ynew = Pnew(:, 2);

  %     trisurf(meshXYpos.triangles, Xnew, Ynew, real(sol.Usol(:, idX)));
  %     hold on;
  %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %   end
  % end
  %
  %
  % for idY = 1:2*opts.numCellsInfinite
  %   sol = load([nomdossier, 'solution_neg_Y_', num2str(idY)]);

  %   for idX = 1:numCellsSemiInfinite_neg
  %     X = meshXYneg.points(:, 1) - (idX - 1) * opts.period;
  %     Y = meshXYneg.points(:, 2) + (idY - opts.numCellsInfinite - 1) * opts.period / opts.cutvec(1);
  %     Pnew = (invT * [X, Y]')';
  %     Xnew = Pnew(:, 1);
  %     Ynew = Pnew(:, 2);

  %     trisurf(meshXYneg.triangles, Xnew, Ynew, real(sol.Usol(:, idX)));
  %     hold on;
  %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %   end
  % end

  % xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
  % ylim([-opts.numCellsInfinite/opts.cutvec(1), opts.numCellsInfinite/opts.cutvec(1)]);
  xlim([-BBx, BBx]);
  ylim([-BBy, BBy]);
  % caxis([-0.25, 0.25]);
  savefig(mafig, [nomdossier, 'fig_solution_2D_', suffix], 'compact');
end

%% Plot 3D solution
if (compute_3D_sol)
  opts.compute3Dsolution = compute_3D_sol;
  numCellsSemiInfinite_pos = ceil(2/opts.period);
  numCellsSemiInfinite_neg = ceil(2/opts.period);
  opts.numCellsInfinite = 2;
  
  opts.suffix = 'pos';
  construct_solution_3D_par_2(true, meshXYpos, meshLineZ, BCstruct.pos, opts);

  opts.suffix = 'neg';
  construct_solution_3D_par_2(true, meshXYneg, meshLineZ, BCstruct.neg, opts);

  %% X constant
  opts.sol3D.points = [ones(meshYZ.numPoints, 1), meshYZ.points];
  opts.sol3D.idCellsX = numCellsSemiInfinite_pos-1;
  opts.sol3D.idCellsY = (-opts.numCellsInfinite:opts.numCellsInfinite-1);
  opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;

  opts.suffix = 'pos';
  construct_solution_3D_par_2(false, meshXYpos, meshLineZ, BCstruct.pos, opts);

  for idX = 1:length(opts.sol3D.idCellsX)     
    mafig3D_X = figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    
    for idY = 1:length(opts.sol3D.idCellsY)
      sol = load([nomdossier, 'solution_3D_pos_X_', num2str(idX), '_Y_', num2str(idY)]);

      Y = opts.sol3D.points(:, 2) + opts.sol3D.idCellsY(idY);
      Z = opts.sol3D.points(:, 3);

      trisurf(meshYZ.triangles, Y, Z, real(sol.Usol));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
    end
    
    savefig(mafig3D_X, [nomdossier, 'fig_solution_3D_pos_X_', int2str(idX)]);
    close(mafig3D_X);
  end

  %% Y constant
  % ////////// %
  % X positive %
  % ////////// %
  meshZX3Dpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesZ, numNodesX);
  opts.sol3D.points = [meshZX3Dpos.points(:, 2), zeros(meshZX3Dpos.numPoints, 1), meshZX3Dpos.points(:, 1)];
  opts.sol3D.idCellsX = (0:numCellsSemiInfinite_pos-1);
  opts.sol3D.idCellsY = -opts.numCellsInfinite;
  opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
  opts.suffix = 'pos';
  construct_solution_3D_par_2(false, meshXYpos, meshLineZ, BCstruct.pos, opts);
  
  mafig3D_Y = cell(length(opts.sol3D.idCellsY));
  for idY = 1:length(opts.sol3D.idCellsY)
    mafig3D_Y{idY} = figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    
    for idX = 1:length(opts.sol3D.idCellsX)
      sol = load([nomdossier, 'solution_3D_pos_X_', num2str(idX), '_Y_', num2str(idY)]);

      X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
      Z = opts.sol3D.points(:, 3);

      trisurf(meshYZ.triangles, Z, X, real(sol.Usol));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
    end
  end
    
  % ////////// %
  % X negative %
  % ////////// %
  meshZX3Dneg = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0 -opts.period], numNodesZ, numNodesX);
  opts.sol3D.points = [meshZX3Dneg.points(:, 2), zeros(meshZX3Dneg.numPoints, 1), meshZX3Dneg.points(:, 1)];
  opts.sol3D.idCellsX = (-numCellsSemiInfinite_neg+1:0);
  opts.sol3D.idCellsY = -opts.numCellsInfinite;
  opts.suffix = 'neg';
  opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
  construct_solution_3D_par_2(false, meshXYneg, meshLineZ, BCstruct.neg, opts);

  for idY = 1:length(opts.sol3D.idCellsY)
    figure(mafig3D_Y{idY});

    for idX = 1:length(opts.sol3D.idCellsX)
      sol = load([nomdossier, 'solution_3D_neg_X_', num2str(idX), '_Y_', num2str(idY)]);

      X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
      Z = opts.sol3D.points(:, 3);

      trisurf(meshYZ.triangles, Z, X, real(sol.Usol));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
    end
    
    savefig(mafig3D_Y{idY}, [nomdossier, 'fig_solution_3D_Y_', int2str(idY)]);
    close(mafig3D_Y{idY});
  end

  
  %% Z constant = 1, X positive
  mafig3D_Z = figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');
  
  % ////////// %
  % X positive %
  % ////////// %
  meshXY3Dpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesX, numNodesY);
  opts.sol3D.points = [meshXY3Dpos.points(:, 1), meshXY3Dpos.points(:, 2), ones(meshXY3Dpos.numPoints, 1)];
  opts.sol3D.idCellsX = (0:numCellsSemiInfinite_pos-1);
  opts.sol3D.idCellsY = (-opts.numCellsInfinite:opts.numCellsInfinite-1);
  opts.suffix = 'pos';
  opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
  construct_solution_3D_par_2(false, meshXYpos, meshLineZ, BCstruct.pos, opts);

  for idY = 1:length(opts.sol3D.idCellsY)
    for idX = 1:length(opts.sol3D.idCellsX)
      sol = load([nomdossier, 'solution_3D_pos_X_', num2str(idX), '_Y_', num2str(idY)]);

      X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
      Y = opts.sol3D.points(:, 2) + opts.sol3D.idCellsY(idY);

      trisurf(meshYZ.triangles, X, Y, real(sol.Usol));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
    end
  end

  % ////////// %
  % X negative %
  % ////////// %
  meshXY3Dneg = FEPack.meshes.MeshRectangle(struct_mesh, [0 -opts.period], [0  opts.period], numNodesX, numNodesY);
  opts.sol3D.points = [meshXY3Dneg.points(:, 1), meshXY3Dneg.points(:, 2), ones(meshXY3Dneg.numPoints, 1)];
  opts.sol3D.idCellsX = (-numCellsSemiInfinite_neg+1:0);
  opts.sol3D.idCellsY = (-opts.numCellsInfinite:opts.numCellsInfinite-1);
  opts.suffix = 'neg';
  opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
  construct_solution_3D_par_2(false, meshXYneg, meshLineZ, BCstruct.neg, opts);

  for idY = 1:length(opts.sol3D.idCellsY)
    for idX = 1:length(opts.sol3D.idCellsX)
      sol = load([nomdossier, 'solution_3D_neg_X_', num2str(idX), '_Y_', num2str(idY)]);

      X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
      Y = opts.sol3D.points(:, 2) + opts.sol3D.idCellsY(idY);

      trisurf(meshYZ.triangles, X, Y, real(sol.Usol));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
    end
  end
  savefig(mafig3D_Z, [nomdossier, 'fig_solution_3D_Z']);
  close(mafig3D_Z);

  system(['rm ', nomdossier, 'solution_3D_*']);
end

% if (compute_3D_sol)
%   opts.compute3Dsolution = compute_3D_sol;
  
%   %% X constant
%   opts.sol3D.points = [ones(meshYZ.numPoints, 1), meshYZ.points];
%   opts.sol3D.idCellsX = numCellsSemiInfinite_pos-1;
%   opts.sol3D.idCellsY = (-opts.numCellsInfinite:opts.numCellsInfinite-1);
%   opts.suffix = 'pos';
%   opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;

%   construct_solution_3D_par(+1, meshXYpos, meshLineZ, BCstruct.pos, opts);

%   for idX = 1:length(opts.sol3D.idCellsX)     
%     mafig3D_X = figure;
%     set(groot,'defaultAxesTickLabelInterpreter','latex');
%     set(groot,'defaulttextinterpreter','latex');
%     set(groot,'defaultLegendInterpreter','latex');
    
%     for idY = 1:length(opts.sol3D.idCellsY)
%       sol = load([nomdossier, 'solution_3D_pos_X_', num2str(idX), '_Y_', num2str(idY)]);

%       Y = opts.sol3D.points(:, 2) + opts.sol3D.idCellsY(idY);
%       Z = opts.sol3D.points(:, 3);

%       trisurf(meshYZ.triangles, Y, Z, real(sol.Usol));
%       hold on;
%       view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%       set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
%     end
    
%     savefig(mafig3D_X, [nomdossier, 'fig_solution_3D_pos_X_', int2str(idX)]);
%     close(mafig3D_X);
%   end

%   %% Y constant
%   opts.sol3D.idCellsY = -opts.numCellsInfinite;
%   for idY = 1:length(opts.sol3D.idCellsY)     
%     % ////////// %
%     % X positive %
%     % ////////// %
%     meshZX3Dpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesZ, numNodesX);
%     opts.sol3D.points = [meshZX3Dpos.points(:, 2), zeros(meshZX3Dpos.numPoints, 1), meshZX3Dpos.points(:, 1)];
%     opts.sol3D.idCellsX = (0:numCellsSemiInfinite_pos-1);
%     opts.suffix = 'pos';
%     opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
%     construct_solution_3D_par(+1, meshXYpos, meshLineZ, BCstruct.pos, opts);

%     mafig3D_Y = figure;
%     set(groot,'defaultAxesTickLabelInterpreter','latex');
%     set(groot,'defaulttextinterpreter','latex');
%     set(groot,'defaultLegendInterpreter','latex');
    
%     for idX = 1:length(opts.sol3D.idCellsX)
%       sol = load([nomdossier, 'solution_3D_pos_X_', num2str(idX), '_Y_', num2str(idY)]);

%       X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
%       Z = opts.sol3D.points(:, 3);

%       trisurf(meshYZ.triangles, Z, X, real(sol.Usol));
%       hold on;
%       view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%       set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
%     end
    
%     % ////////// %
%     % X negative %
%     % ////////// %
%     meshZX3Dneg = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0 -opts.period], numNodesZ, numNodesX);
%     opts.sol3D.points = [meshZX3Dneg.points(:, 2), zeros(meshZX3Dneg.numPoints, 1), meshZX3Dneg.points(:, 1)];
%     opts.sol3D.idCellsX = (-numCellsSemiInfinite_neg+1:0);
%     opts.suffix = 'neg';
%     opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
%     construct_solution_3D_par(+1, meshXYneg, meshLineZ, BCstruct.neg, opts);

%     for idX = 1:length(opts.sol3D.idCellsX)
%       sol = load([nomdossier, 'solution_3D_neg_X_', num2str(idX), '_Y_', num2str(idY)]);

%       X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
%       Z = opts.sol3D.points(:, 3);

%       trisurf(meshYZ.triangles, Z, X, real(sol.Usol));
%       hold on;
%       view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%       set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
%     end
    
%     savefig(mafig3D_Y, [nomdossier, 'fig_solution_3D_Y_', int2str(idY)]);
%     close(mafig3D_Y);
%   end

  
%   %% Z constant = 1, X positive
%   mafig3D_Z = figure;
%   set(groot,'defaultAxesTickLabelInterpreter','latex');
%   set(groot,'defaulttextinterpreter','latex');
%   set(groot,'defaultLegendInterpreter','latex');
  
%   % ////////// %
%   % X positive %
%   % ////////// %
%   meshXY3Dpos = FEPack.meshes.MeshRectangle(struct_mesh, [0  opts.period], [0  opts.period], numNodesX, numNodesY);
%   opts.sol3D.points = [meshXY3Dpos.points(:, 1), meshXY3Dpos.points(:, 2), ones(meshXY3Dpos.numPoints, 1)];
%   opts.sol3D.idCellsX = (0:numCellsSemiInfinite_pos-1);
%   opts.sol3D.idCellsY = (-opts.numCellsInfinite:opts.numCellsInfinite-1);
%   opts.suffix = 'pos';
%   opts.numCellsSemiInfinite = numCellsSemiInfinite_pos;
%   construct_solution_3D_par(+1, meshXYpos, meshLineZ, BCstruct.pos, opts);

%   for idY = 1:length(opts.sol3D.idCellsY)
%     for idX = 1:length(opts.sol3D.idCellsX)
%       sol = load([nomdossier, 'solution_3D_pos_X_', num2str(idX), '_Y_', num2str(idY)]);

%       X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
%       Y = opts.sol3D.points(:, 2) + opts.sol3D.idCellsY(idY);

%       trisurf(meshYZ.triangles, X, Y, real(sol.Usol));
%       hold on;
%       view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%       set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
%     end
%   end

%   % ////////// %
%   % X negative %
%   % ////////// %
%   meshXY3Dneg = FEPack.meshes.MeshRectangle(struct_mesh, [0 -opts.period], [0  opts.period], numNodesX, numNodesY);
%   opts.sol3D.points = [meshXY3Dneg.points(:, 1), meshXY3Dneg.points(:, 2), ones(meshXY3Dneg.numPoints, 1)];
%   opts.sol3D.idCellsX = (-numCellsSemiInfinite_neg+1:0);
%   opts.sol3D.idCellsY = (-opts.numCellsInfinite:opts.numCellsInfinite-1);
%   opts.suffix = 'neg';
%   opts.numCellsSemiInfinite = numCellsSemiInfinite_neg;
%   construct_solution_3D_par(+1, meshXYneg, meshLineZ, BCstruct.neg, opts);

%   for idY = 1:length(opts.sol3D.idCellsY)
%     for idX = 1:length(opts.sol3D.idCellsX)
%       sol = load([nomdossier, 'solution_3D_neg_X_', num2str(idX), '_Y_', num2str(idY)]);

%       X = opts.sol3D.points(:, 1) + opts.sol3D.idCellsX(idX);
%       Y = opts.sol3D.points(:, 2) + opts.sol3D.idCellsY(idY);

%       trisurf(meshYZ.triangles, X, Y, real(sol.Usol));
%       hold on;
%       view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%       set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
%     end
%   end
%   savefig(mafig3D_Z, [nomdossier, 'fig_solution_3D_Z']);
%   close(mafig3D_Z);

%   system(['rm ', nomdossier, 'solution_3D_*']);
% end

%%
system(['rm ', nomdossier, 'half_guide_*']);
system(['rm ', nomdossier, 'local_cell_sol_*']);
system(['rm ', nomdossier, 'sol_trace_Floquet_*']);
if (delete_sol)
  system(['rm ', nomdossier, 'solution_*']);
end
