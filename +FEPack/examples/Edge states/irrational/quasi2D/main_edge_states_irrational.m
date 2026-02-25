%% Import FEPack
clear; clc;
import FEPack.*

%% Problem-related parameters
HcObj = tools.HoneycombObject('none', 'none');

% Honeycomb lattice potentials
% centers = [-1/sqrt(3), 1/sqrt(3); 0, 0];
% ampsV   = [0;  0];
% ampsW   = [10; 10];
% rads    = [0.2; 0.2];
% HcObj.V = @(x) FEPack.tools.FUN_atomic_potential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsV, rads);
% HcObj.W = @(x) 1 + FEPack.tools.FUN_atomic_potential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsW, rads);

HcObj.V = @(x) 10 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
               10 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
               10 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

HcObj.W = @(x)      sin(x(:, 1:2) *  HcObj.dualVec1) +...
                    sin(x(:, 1:2) *  HcObj.dualVec2) +...
                    sin(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

% Edge
HcObj.edge.a1 =  1;
HcObj.edge.b1 =  sqrt(2);

% Domain wall 
HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);
pbinputs.delta = 10;

%% Energy/parallel quasi-momentum range
pbinputs.minKpar = -pi;
pbinputs.maxKpar = +pi;
pbinputs.numKpar = 16;
pbinputs.center_parallel_quasi_momentum_around_K = false;

% Energy range
pbinputs.minEgy = 11.7341-10; % 17.5460-1e-1;
pbinputs.maxEgy = 11.7341+10;% 17.5460+1e-1;
pbinputs.numEgy = 16;

%% Meshes and boundary conditions
% # nodes per X/Y edge for periodicity cells
pbinputs.numNodesX      =  8;
pbinputs.numNodesY      =  8; 
pbinputs.numNodesZ      =  8;
pbinputs.numNodesXint   =  8;
pbinputs.numNodesYint   =  8;
pbinputs.quadratureRatio = 8;
pbinputs.is_mesh_structured = 1;

% Basis functions for boundary conditions ("Lagrange" for now)
pbinputs.bc_basis_functions_X = "Lagrange";
pbinputs.bc_basis_functions_S = "Lagrange";

%% Misc. 
% Outputs folder
opts.folder_name = 'outputs/';

% Should you use parallel computing?
opts.parallel_use = false;
opts.parallel_numcores = 4;

% plot the coefficients
opts.plot_coefficients = false;

% Plot dispersion functions, save them
opts.plot_disp = false;
opts.save_disp = true;

%% Initialization

% Construct operator structures
[op_pos, op_neg, cutvec] = construct_ops(pbinputs, HcObj);

% Construct FE matrices and boundary-related basis structures
[meshXYpos, meshXYneg, meshYZ, meshLineZ, init_BCstruct_pos, init_BCstruct_neg]...
  = initialize_edge_states_irrational(...
      pbinputs.numNodesX, pbinputs.numNodesY, pbinputs.numNodesZ, pbinputs.is_mesh_structured,...
      pbinputs.bc_basis_functions_X, pbinputs.bc_basis_functions_S,...
      op_pos, op_neg, cutvec,...
      opts.folder_name,...
      opts.parallel_use, opts.parallel_numcores...
    );

%% Launch it!
numKcell = 1;
numEcell = 1;

kpars = linspace(pbinputs.minKpar, pbinputs.maxKpar, numKcell+1);
Egies = linspace(pbinputs.minEgy,  pbinputs.maxEgy,  numEcell+1);

Kstep = kpars(2) - kpars(1);
Estep = Egies(2) - Egies(1);

if (opts.plot_disp)
  figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');
  
  hold on;
end


for idK = 1:numKcell
  for idE = 1:numEcell

    clc;
    pbinputs.minKpar = kpars(idK);
    pbinputs.maxKpar = kpars(idK+1);
    pbinputs.minEgy  = Egies(idE);
    pbinputs.maxEgy  = Egies(idE+1);

    BBkpar = [pbinputs.minKpar, pbinputs.maxKpar];
    BBlambda = [pbinputs.minEgy, pbinputs.maxEgy];

    %% Computations
    [dispmesh, val] = edge_states_irrational(...
      BBkpar, pbinputs.numKpar,...
      BBlambda, pbinputs.numEgy,...
      meshXYpos, meshXYneg, meshYZ, meshLineZ, pbinputs.quadratureRatio,...
      cutvec,...
      init_BCstruct_pos, init_BCstruct_neg,...
      opts.folder_name,...
      opts.parallel_use, opts.parallel_numcores);

    %% Plot the result
    if (opts.plot_disp)
      trisurf(dispmesh.triangles, dispmesh.points(:, 1), dispmesh.points(:, 2), log(abs(val))); hold on;
      % plot3( HcObj.highSymK'*HcObj.vecPer1*[1 1], [pbinputs.minEgy, pbinputs.maxEgy], [1e8 1e8], 'w--');
      % plot3(-HcObj.highSymK'*HcObj.vecPer1*[1 1], [pbinputs.minEgy, pbinputs.maxEgy], [1e8 1e8], 'w--');
    end

    %% Save the result
    if (opts.save_disp)
      save(['outputs/disp_', num2str(idK), '_', num2str(idE)], 'dispmesh', 'val', '-v7.3');
    end
    
  end
end

if (opts.plot_disp)
  axis([kpars(1) kpars(end) Egies(1) Egies(end)]);
  shading interp;
  colormap jet;
  view(2);
  set(gca, 'FontSize', 16);
  colorbar('TickLabelInterpreter', 'latex');
  xlabel('$k_\parallel$');
  ylabel('$E$');
end


%% Auxiliary functions
function [op_pos, op_neg, cutvec] = construct_ops(pbinputs, HcObj)

  % Make sure the edge coefficients are coprime integers
  a1 = HcObj.edge.a1;
  b1 = HcObj.edge.b1;
  a2 = 0;
  b2 = 1;

  % Edge vectors
  edge_vec1      =  a1 * HcObj.vecPer1  + b1 * HcObj.vecPer2;
  edge_vec2      =  a2 * HcObj.vecPer1  + b2 * HcObj.vecPer2;
  edge_dual_vec1 =  b2 * HcObj.dualVec1 - a2 * HcObj.dualVec2;
  edge_dual_vec2 = -b1 * HcObj.dualVec1 + a1 * HcObj.dualVec2;

  % Cut vector
  cutvec = [1; HcObj.edge.b1];

  % Transformation matrices
  Rmat = [[edge_vec2, edge_vec1, zeros(2, 1)]; [0, b1 1]];
  Tmat = [edge_dual_vec2'; edge_dual_vec1'] / (2*pi); 
  Tvec = Tmat' * [1; 0];

  % Potentials
  Vpot3Dinit = @(x) HcObj.V(x(:, 1:2) - x(:, 3) * edge_vec2');
  Wpot3Dinit = @(x) HcObj.W(x(:, 1:2) - x(:, 3) * edge_vec2');

  Vpot3D = @(x) Vpot3Dinit((Rmat * x')');
  Wpot3D = @(x) Wpot3Dinit((Rmat * x')');

  % (+) half-guide
  op_pos.Tmat   = Tmat;
  op_pos.vect   = Tvec;
  op_pos.funQ3D = @(x) Vpot3D(x) + pbinputs.delta * Wpot3D(x);
  op_pos.funP3D = 1;
  op_pos.funR3D = 1;

  % (-) half-guide
  op_neg.Tmat   = Tmat;
  op_neg.vect   = Tvec;
  op_neg.funQ3D = @(x) Vpot3D(x) - pbinputs.delta * Wpot3D(x);
  op_neg.funP3D = 1;
  op_neg.funR3D = 1;

  % % Interior domain
  % yInt = ceil(1/(2*pi*pbinputs.delta));
  % BBint = [-yInt, yInt];
  % op_int.Tmat = Tmat;
  % op_int.vect = Tvec;
  % op_int.funQ = @(x) Vpot(x) + pbinputs.delta * HcObj.kappa(2*pi*pbinputs.delta * x(:, 1)) .* Wpot(x);
  % op_int.funP = 1;
  % op_int.funR = 1;
end
