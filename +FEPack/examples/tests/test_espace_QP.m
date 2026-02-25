clear; clc;
import FEPack.*

% profile ON
opts.omega = 5 + 0.25i;% 5 + 0.25i;
cutangle = pi/4;% pi/3;
cutslope = tan(cutangle);
matN = [1, -cutslope; 0, 1];
opts.cutmat = [1, 0; 0, 1; 0, cutslope] * matN;

funumNodes3D = @(fun2D, x) fun2D([x(:, 1) + x(:, 3), x(:, 2), zeros(size(x, 1), 1)]);

mu2Dpos = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho2Dpos = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
mu3Dpos = @(x) funumNodes3D(mu2Dpos, x);
rho3Dpos = @(x) funumNodes3D(rho2Dpos, x);

mu2Dneg = @(x) ones(size(x, 1), 1); % @(x) 1 + 0.5*sin(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
rho2Dneg = @(x) ones(size(x, 1), 1); % @(x) 1 + 0.25*sin(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
mu3Dneg = @(x) funumNodes3D(mu2Dneg, x);
rho3Dneg = @(x) funumNodes3D(rho2Dneg, x);

mu3Dint =  @(x) (x(:, 1) >= 0) .*  mu3Dpos(x) + (x(:, 1) <  0) .* mu3Dneg(x);
rho3Dint = @(x) (x(:, 1) >= 0) .* rho3Dpos(x) + (x(:, 1) <  0) .* rho3Dneg(x);

O = [-0.5, 0.5];
F2Ddroit = @(x) 100*FEPack.tools.cutoff(sqrt((x(:, 1)-O(1)).^2 + (x(:, 2)-O(2)).^2), -0.3, 0.3);
% F2Ddroit = @(x) FEPack.tools.cutoff(x(:, 1)-O(1), -0.3, 0.3) .* FEPack.tools.cutoff(x(:, 2)-O(2), -0.3, 0.3);
% F2D = @(x) F2Ddroit((matN * x')');
F3D = @(x) F2Ddroit(x);

% Supp
intBounds.pos = +0.5;
intBounds.neg = -1;
structmesh = 0;
basis_functions = 'Fourier';
u = pdes.PDEObject; v = dual(u);
% volBilinearIntg = @(muco, rhoco) (muco * grad(u)) * grad(v) - (opts.omega^2) * ((rhoco*id(u))*id(v));
volBilinearIntg = @(muco, rhoco) (muco * gradDir(u, opts.cutmat)) * gradDir(v, opts.cutmat) - (opts.omega^2) * ((rhoco*id(u))*id(v));
plot_coefficients = false;
straight_boundary = true;
compareU = true;

numNodes2D = 32;
numNodes3D = 16;

%% Parameters for the positive half-guide
%  //////////////////////////////////////
mesh2Dpos = meshes.MeshRectangle(structmesh, [0 1], [0 1], numNodes2D, numNodes2D);
mesh3Dpos = meshes.MeshCuboid(structmesh, [0 1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmin'));
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmax'));
else
  FourierIds = [0, floor(numNodes3D/4), 0];
  BCstruct_pos.spB0 = spaces.FourierBasis(mesh3Dpos.domain('xmin'), FourierIds);
  BCstruct_pos.spB1 = spaces.FourierBasis(mesh3Dpos.domain('xmax'), FourierIds);
end

BCstruct_pos.BCdu = 0.0;
BCstruct_pos.BCu = 1.0;% 1.0; % 1i*opts.omega;% 1.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_pos.representation = '';
volBilinearIntg_pos = volBilinearIntg(@(x)  mu3Dpos(x + sparse(1:size(x, 1), 1, intBounds.pos, size(x, 1), size(x, 2))),...
                                      @(x) rho3Dpos(x + sparse(1:size(x, 1), 1, intBounds.pos, size(x, 1), size(x, 2))));

%% Parameters for the negative half-guide
%  //////////////////////////////////////
mesh2Dneg = meshes.MeshRectangle(structmesh, [0 -1], [0 1], numNodes2D, numNodes2D);
mesh3Dneg = meshes.MeshCuboid(structmesh, [0 -1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmin'));
  BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmax'));
else
  FourierIds = [0, floor(numNodes3D/4), 0];
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh3Dneg.domain('xmin'), FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh3Dneg.domain('xmax'), FourierIds);
end

BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;% 0.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_neg.representation = '';
numCells_neg = 4;
volBilinearIntg_neg = volBilinearIntg(@(x)  mu3Dneg(x + sparse(1:size(x, 1), 1, intBounds.neg, size(x, 1), size(x, 2))),...
                                      @(x) rho3Dneg(x + sparse(1:size(x, 1), 1, intBounds.neg, size(x, 1), size(x, 2))));

%% Parameters for the interior problem
%  //////////////////////////////////////
Lint = abs(intBounds.pos - intBounds.neg);
mesh2Dint = meshes.MeshRectangle(structmesh, [intBounds.neg, intBounds.pos], [0 1], floor(Lint*numNodes2D), numNodes2D);
mesh3Dint = meshes.MeshCuboid(structmesh, [intBounds.neg, intBounds.pos], [0 1], [0 1], floor(Lint*numNodes3D), numNodes3D, numNodes3D);

if strcmpi(basis_functions, 'Lagrange')
  spBint_neg = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dint.domain('xmin'));
  spBint_pos = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dint.domain('xmax'));
else
  FourierIds = [0, floor(numNodes3D/4), 0];
  spBint_pos = spaces.FourierBasis(mesh3Dint.domain('xmin'), FourierIds);
  spBint_neg = spaces.FourierBasis(mesh3Dint.domain('xmax'), FourierIds);
end

volBilinearIntg_int = volBilinearIntg(mu3Dint, rho3Dint);
volLinearIntg = F3D * id(v);

numCellsSemiInfinite_pos = 5;
numCellsSemiInfinite_neg = 5;
numCellsInfinite = 4;
numFloquetPoints = 20;

%% Plot the coefficients and the source term
if (plot_coefficients)
  set(groot,'defaultAxesTickLabelInterpreter','latex'); %#ok
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  matNinv = [1, cutslope; 0, 1];
  mu2Dint =  @(x) (x(:, 1) >= 0) .*  mu2Dpos((matNinv * x')') + (x(:, 1) <  0) .* mu2Dneg((matNinv * x')');
  rho2Dint = @(x) (x(:, 1) >= 0) .* rho2Dpos((matNinv * x')') + (x(:, 1) <  0) .* rho2Dneg((matNinv * x')');

  for idS = 1:(2*numCellsSemiInfinite_pos)
    for idI = 1:(2*numCellsInfinite)
      X = mesh2Dpos.points(:, 1) + (idS - numCellsSemiInfinite_pos - 1);
      Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);
      figure(1);
      trisurf(mesh2Dpos.triangles, X, Y, mu2Dint([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

      figure(2);
      trisurf(mesh2Dpos.triangles, X, Y, rho2Dint([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

      figure(3);
      trisurf(mesh2Dpos.triangles, X, Y, F2Ddroit([X, Y]));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end
end

%%
% Compute guide solution
U3D = PeriodicSpaceBVP(1, 2,...
                       volBilinearIntg_pos, mesh3Dpos, BCstruct_pos, numCellsSemiInfinite_pos,...
                       volBilinearIntg_neg, mesh3Dneg, BCstruct_neg, numCellsSemiInfinite_neg,...
                       volBilinearIntg_int, volLinearIntg, mesh3Dint, spBint_pos, spBint_neg,...
                       numCellsInfinite, numFloquetPoints, opts);

%%
% Take the trace
% Positive side
% /////////////
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
%%
% % Negative side
% % /////////////
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

%%
% Interior domain
% ///////////////
N2Dint = mesh2Dint.numPoints;
U2D.interior = zeros(N2Dint, size(U3D.interior, 2));
dom = mesh3Dint.domain('volumic');
for idI = 1:2*numCellsInfinite
  Icell = idI;
  X = mesh2Dint.points(:, 1);
  Y = mesh2Dint.points(:, 2); % ones(mesh2Dint.numPoints, 1);
  Z = FEPack.tools.mymod(cutslope * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.interior(:, Icell) = reshape(sum(reshape(coos .* U3D.interior(elts, Icell), dom.dimension+1, []), 1), N2Dint, []);
end


%% Plot U
if straight_boundary
  Tplot = eye(2);
else
  Tplot = [1, cutslope; 0, 1];
end

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
if (compareU)
  subplot(1, 2, 1);
end
for idS = 1:numCellsSemiInfinite_pos
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
    X = mesh2Dpos.points(:, 1) + (idS - 1) + intBounds.pos;
    Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);

    trisurf(mesh2Dpos.triangles, Tplot(1,1)*X + Tplot(1,2)*Y, Tplot(2,1)*X + Tplot(2,2)*Y, real(U2D.positive(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
  end
end
%%
for idS = 1:numCellsSemiInfinite_neg
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
    X = mesh2Dneg.points(:, 1) - (idS - 1) + intBounds.neg;
    Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
    trisurf(mesh2Dneg.triangles, Tplot(1,1)*X + Tplot(1,2)*Y, Tplot(2,1)*X + Tplot(2,2)*Y, real(U2D.negative(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
  end
end
%%
% figure;
for idI = 1:size(U2D.interior, 2)%2*numCellsInfinite
  Icell = idI;
  X = mesh2Dint.points(:, 1);
  Y = mesh2Dint.points(:, 2) + (idI - numCellsInfinite - 1);
  trisurf(mesh2Dint.triangles, Tplot(1,1)*X + Tplot(1,2)*Y, Tplot(2,1)*X + Tplot(2,2)*Y, real(U2D.interior(:, Icell)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  % colormap jet
end

xlim([intBounds.neg - numCellsSemiInfinite_neg + 1, intBounds.pos + numCellsSemiInfinite_pos - 1]);
ylim([-numCellsInfinite, numCellsInfinite]);

%% Compare U in the rational case
if (compareU)

  matNinv = [1, cutslope; 0, 1];

  mu2DposD = @(x) mu2Dpos([x(:, 1) + cutslope*x(:, 2), x(:, 2)]);
  rho2DposD = @(x) rho2Dpos([x(:, 1) + cutslope*x(:, 2), x(:, 2)]);
  BCstruct2Dpos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmin'));
  BCstruct2Dpos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmax'));
  BCstruct2Dpos.BCdu = BCstruct_pos.BCdu;
  BCstruct2Dpos.BCu = BCstruct_pos.BCu;
  volBilinearIntg2Dpos = volBilinearIntg(@(x)  mu2DposD(x + sparse(1:size(x, 1), 1, intBounds.pos, size(x, 1), size(x, 2))),...
                                         @(x) rho2DposD(x + sparse(1:size(x, 1), 1, intBounds.pos, size(x, 1), size(x, 2))));

  %
  mu2DnegD = @(x) mu2Dneg([x(:, 1) + cutslope*x(:, 2), x(:, 2)]);
  rho2DnegD = @(x) rho2Dneg([x(:, 1) + cutslope*x(:, 2), x(:, 2)]);
  BCstruct2Dneg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmin'));
  BCstruct2Dneg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmax'));
  BCstruct2Dneg.BCdu = BCstruct_neg.BCdu;
  BCstruct2Dneg.BCu = BCstruct_neg.BCu;
  volBilinearIntg2Dneg = volBilinearIntg(@(x)  mu2DnegD(x + sparse(1:size(x, 1), 1, intBounds.neg, size(x, 1), size(x, 2))),...
                                         @(x) rho2DnegD(x + sparse(1:size(x, 1), 1, intBounds.neg, size(x, 1), size(x, 2))));

  %
  mu2DintD =  @(x) (x(:, 1) >= 0) .*  mu2DposD(x) + (x(:, 1) <  0) .*  mu2DnegD(x);
  rho2DintD = @(x) (x(:, 1) >= 0) .* rho2DposD(x) + (x(:, 1) <  0) .* rho2DnegD(x);
  spBint2Dneg = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dint.domain('xmin'));
  spBint2Dpos = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dint.domain('xmax'));
  volBilinearIntg2Dint = volBilinearIntg(mu2DintD, rho2DintD);

  %
  Ur = PeriodicSpaceBVP(1, 2,...
                        volBilinearIntg2Dpos, mesh2Dpos, BCstruct2Dpos, numCellsSemiInfinite_pos,...
                        volBilinearIntg2Dneg, mesh2Dneg, BCstruct2Dneg, numCellsSemiInfinite_neg,...
                        volBilinearIntg2Dint, volLinearIntg, mesh2Dint, spBint2Dpos, spBint2Dneg,...
                        numCellsInfinite, numFloquetPoints, opts);
  %%
  % figure;
  % figure(2)
  subplot(1, 2, 2);
  
  for idS = 1:numCellsSemiInfinite_pos
    for idI = 1:2*numCellsInfinite
      Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
      X = mesh2Dpos.points(:, 1) + (idS - 1) + intBounds.pos;
      Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);

      trisurf(mesh2Dpos.triangles, X, Y, real(Ur.positive(:, Icell)));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      % colormap jet
    end
  end
  %
  for idS = 1:numCellsSemiInfinite_neg
    for idI = 1:2*numCellsInfinite
      Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
      X = mesh2Dneg.points(:, 1) - (idS - 1) + intBounds.neg;
      Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
      trisurf(mesh2Dneg.triangles, X, Y, real(Ur.negative(:, Icell)));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      % colormap jet
    end
  end
  %
  % figure;
  for idI = 1:size(U2D.interior, 2)%2*numCellsInfinite
    Icell = idI;
    X = mesh2Dint.points(:, 1);
    Y = mesh2Dint.points(:, 2) + (idI - numCellsInfinite - 1);
    trisurf(mesh2Dint.triangles, X, Y, real(Ur.interior(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
  end
  %
  xlim([intBounds.neg - numCellsSemiInfinite_neg + 1, intBounds.pos + numCellsSemiInfinite_pos - 1]);
  ylim([-numCellsInfinite, numCellsInfinite]);

end
