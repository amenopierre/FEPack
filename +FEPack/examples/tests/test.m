clear; clc;
%%
import FEPack.*
% profile ON
opts.omega = 8 + 0.1i;% 5 + 0.25i;
vecper = [-sqrt(2), 1];% [-sqrt(2), 1];
cutvec = [-vecper(1)/vecper(2); 1.0/vecper(2)];
opts.cutmat = [[1; 0; 0], [0; cutvec]];
cutslope = cutvec(2) / cutvec(1);


mu2DposCart = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho2DposCart = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
mu2Dpos = @(x) mu2DposCart([x(:, 1) + cutvec(1)*x(:, 2), cutvec(2)*x(:, 2), zeros(size(x, 1), 1)]);
rho2Dpos = @(x) rho2DposCart([x(:, 1) + cutvec(1)*x(:, 2), cutvec(2)*x(:, 2), zeros(size(x, 1), 1)]);

mu2DnegCart = @(x) ones(size(x, 1), 1); % @(x) 1 + 0.5*sin(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
rho2DnegCart = @(x) ones(size(x, 1), 1); % @(x) 1 + 0.25*sin(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
mu2Dneg = @(x) mu2DnegCart([x(:, 1) + cutvec(1)*x(:, 2), cutvec(2)*x(:, 2), zeros(size(x, 1), 1)]);
rho2Dneg = @(x) rho2DnegCart([x(:, 1) + cutvec(1)*x(:, 2), cutvec(2)*x(:, 2), zeros(size(x, 1), 1)]);

mu2Dint =  @(x) (x(:, 1) >= 0) .*  mu2Dpos(x) + (x(:, 1) <  0) .* mu2Dneg(x);
rho2Dint = @(x) (x(:, 1) >= 0) .* rho2Dpos(x) + (x(:, 1) <  0) .* rho2Dneg(x);

G = @(x) FEPack.tools.cutoff(x(:, 2), -0.3, 0.3);


% Supp
structmesh = 0;
basis_functions = 'Fourier';
u = pdes.PDEObject; v = dual(u);
% volBilinearIntg = @(muco, rhoco) (muco * grad(u)) * grad(v) - (opts.omega^2) * ((rhoco*id(u))*id(v));
volBilinearIntg = @(muco, rhoco) (muco * gradDir(u, opts.cutmat)) * gradDir(v, opts.cutmat) - (opts.omega^2) * ((rhoco*id(u))*id(v));

volBilinearIntg2D = volBilinearIntg(mu2Dpos, rho2Dpos);
jumpLinearIntg = G * id(v);

plot_coefficients = false;
compareU = false;

numNodes2D = 64;

mesh2D = meshes.MeshRectangle(structmesh, [-13 13], [-13 13], numNodes2D, numNodes2D);

AA = FEPack.pdes.Form.intg(





%% Parameters for the positive half-guide
%  //////////////////////////////////////
mesh3Dpos = meshes.MeshCuboid(structmesh, [0 1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmin'));
  BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmax'));
else
  FourierIds = [0, floor(numNodes3D/2), 0];
  BCstruct_pos.spB0 = spaces.FourierBasis(mesh3Dpos.domain('xmin'), FourierIds);
  BCstruct_pos.spB1 = spaces.FourierBasis(mesh3Dpos.domain('xmax'), FourierIds);
end

BCstruct_pos.BCdu = 0.0;
BCstruct_pos.BCu = 1.0;% 1.0; % 1i*opts.omega;% 1.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_pos.representation = 'weak evaluation';
volBilinearIntg_pos = volBilinearIntg(mu3Dpos, rho3Dpos);

%% Parameters for the negative half-guide
%  //////////////////////////////////////
mesh2Dneg = meshes.MeshRectangle(structmesh, [0 -1], [0 1], numNodes2D, numNodes2D);
mesh3Dneg = meshes.MeshCuboid(structmesh, [0 -1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

if strcmpi(basis_functions, 'Lagrange')
  BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmin'));
  BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmax'));
else
  FourierIds = [0, floor(numNodes3D/2), 0];
  BCstruct_neg.spB0 = spaces.FourierBasis(mesh3Dneg.domain('xmin'), FourierIds);
  BCstruct_neg.spB1 = spaces.FourierBasis(mesh3Dneg.domain('xmax'), FourierIds);
end

BCstruct_neg.BCdu = 0.0;
BCstruct_neg.BCu = 1.0;% 0.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct_neg.representation = 'weak evaluation';
numCells_neg = 4;
volBilinearIntg_neg = volBilinearIntg(mu3Dneg, rho3Dneg);

%% Parameters for the interface problem
%  //////////////////////////////////////
jumpLinearIntg = G3D * id(v);

numCellsSemiInfinite_pos = 5;
numCellsSemiInfinite_neg = 5;
numCellsInfinite = 4;
numFloquetPoints = 50;

%% Plot the coefficients and the source term
if (false)
  set(groot,'defaultAxesTickLabelInterpreter','latex'); %#ok
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  mu2Dint =  @(x) (x(:, 1) >= 0) .*  mu2Dpos(x) + (x(:, 1) <  0) .* mu2Dneg(x);
  rho2Dint = @(x) (x(:, 1) >= 0) .* rho2Dpos(x) + (x(:, 1) <  0) .* rho2Dneg(x);

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

    end
  end

  figure(3);
  x = [zeros(256, 1), linspace(-2, 2, 256)', zeros(256, 1)];
  plot(x(:, 2), G(x));
end

%%
% Compute guide solution
U3D = PeriodicSpaceJumpBVP(1, 2,...
                           volBilinearIntg_pos, mesh3Dpos, BCstruct_pos, numCellsSemiInfinite_pos,...
                           volBilinearIntg_neg, mesh3Dneg, BCstruct_neg, numCellsSemiInfinite_neg,...
                           jumpLinearIntg, numCellsInfinite, numFloquetPoints, opts);

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


%% Plot U
figure;%(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
if (compareU)
  subplot(1, 2, 1);
end
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

xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
ylim([-numCellsInfinite, numCellsInfinite]);

% figure('Position', get(0, 'Screensize'), 'visible', 'off');
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off;


%% Compare U in the rational case
if (compareU)
  
  BCstruct2Dpos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmin'));
  BCstruct2Dpos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmax'));
  BCstruct2Dpos.BCdu = BCstruct_pos.BCdu;
  BCstruct2Dpos.BCu = BCstruct_pos.BCu;
  volBilinearIntg2Dpos = volBilinearIntg(mu2Dpos, rho2Dpos);

  %
  BCstruct2Dneg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmin'));
  BCstruct2Dneg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmax'));
  BCstruct2Dneg.BCdu = BCstruct_neg.BCdu;
  BCstruct2Dneg.BCu = BCstruct_neg.BCu;
  volBilinearIntg2Dneg = volBilinearIntg(mu2Dneg, rho2Dneg);

  %
  jumpLinearIntg2D = G * id(v);
  
  %
  Ur = PeriodicSpaceJumpBVP(1, 2,...
                        volBilinearIntg2Dpos, mesh2Dpos, BCstruct2Dpos, numCellsSemiInfinite_pos,...
                        volBilinearIntg2Dneg, mesh2Dneg, BCstruct2Dneg, numCellsSemiInfinite_neg,...
                        jumpLinearIntg2D, numCellsInfinite, numFloquetPoints, opts);
  %%
  % figure;
  % figure(2)
  subplot(1, 2, 2);
  
  for idS = 1:numCellsSemiInfinite_pos
    for idI = 1:2*numCellsInfinite
      Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
      X = mesh2Dpos.points(:, 1) + (idS - 1);
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
      X = mesh2Dneg.points(:, 1) - (idS - 1);
      Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
      trisurf(mesh2Dneg.triangles, X, Y, real(Ur.negative(:, Icell)));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      % colormap jet
    end
  end
  %
  %
  % xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
  % ylim([-numCellsInfinite, numCellsInfinite]);

end

