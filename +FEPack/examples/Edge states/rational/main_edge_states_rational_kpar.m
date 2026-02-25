% function main_edge_states_rational_kpar(a1, b1)
  %% Import FEPack
  import FEPack.*
  % if (nargin < 2), b1 = 0; end
  % if (nargin < 1), a1 = 1; end
  clear; clc;
  a1 = 1; b1 = -1;
  

  %% Problem-related parameters
  HcObj = tools.HoneycombObject('none', 'none', 'none');

  % Honeycomb lattice potentials
  HcObj.V = @(x) 20 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
                 20 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
                 20 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

  HcObj.W = @(x)  0 * sin(x(:, 1:2) *  HcObj.dualVec1) +...
                  0 * sin(x(:, 1:2) *  HcObj.dualVec2) +...
                  0 * sin(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

  HcObj.A = @(x)  0 +...
                  1 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
                  1 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
                  1 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

  % Edge
  HcObj.edge.a1 = a1;
  HcObj.edge.b1 = b1;
  
  % Problem type ('interior' or 'interface')
  pbinputs.problem_type = 'interior';

  % Domain wall 
  HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);
  pbinputs.delta = 0.5;

  %% Energy/parallel quasi-momentum range
  pbinputs.minKpar = -pi;
  pbinputs.maxKpar = +pi;
  pbinputs.center_parallel_quasi_momentum_around_K = false;

  % Energy range
  % Dirac point 11.7341
  pbinputs.minEgy = -8; % 11.7341-10; % 17.5460-1e-1;
  pbinputs.maxEgy = 12; % 11.7341+10; % 17.5460+1e-1;
  
  %% Meshes and boundary conditions
  % # nodes per X/Y edge for periodicity cells
  numNodes = 16;
  pbinputs.numNodesXpos = numNodes;
  pbinputs.numNodesYpos = numNodes; 
  pbinputs.numNodesXneg = numNodes;
  pbinputs.numNodesYneg = numNodes;
  pbinputs.numNodesXint = numNodes;
  pbinputs.numNodesYint = numNodes;

  % Basis functions for boundary conditions ("Lagrange" or "Fourier")
  pbinputs.bc_basis_functions_pos = "Lagrange";
  pbinputs.bc_basis_functions_neg = "Lagrange";

  %% Misc. 
  % Should you use parallel computing?
  parallel_use = false;
  parallel_numcores = 2;

  % plot the coefficients
  plot_coefficients = false;

  %% Initialization
  [mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = construct_mesh_and_FEmatrices(pbinputs, HcObj);

  %% Compute indicator for edge state curves
  numKpoints = 32;
  kpars = linspace(-pi, pi, numKpoints);
  minEgy = pbinputs.minEgy;
  maxEgy = pbinputs.maxEgy;
  problem_type = pbinputs.problem_type;
  
  rootBB = [minEgy, maxEgy];
  maxDepth = 12;
  force_depth_gap  = 8;
  force_depth_band = 6;
  force_depth_gap  = min(force_depth_gap,  maxDepth);
  force_depth_band = min(force_depth_band, maxDepth);
  threshold = 1e4;

  inside_band    = @(fvals) all(isnan(fvals), 2);
  on_band_border = @(fvals) (any(isnan(fvals), 2) & ~inside_band(fvals));
  on_edge_border = @(fvals) (any(abs(fvals) > threshold, 2));
  early_split    = @(depth, force_depth) (depth <= force_depth);

  % The rule: split leaves on the border of continuous spectrum XOR
  % min depth hasn't been reached yet or the indicator function is 
  % above a threshold.
  rules = @(fvals, depth) ( on_band_border(fvals)) |...
                          ( on_edge_border(fvals)) |...
                          ( inside_band(fvals) & early_split(depth, force_depth_band)) |...
                          (~inside_band(fvals) & early_split(depth, force_depth_gap));
  
  for idK = 1:numKpoints
    clc; disp(idK)
    kpar = kpars(idK);

    bt = FEPack.meshes.BinaryTree(maxDepth, rootBB, rules);

    indicator_function = @(P) edge_states_rational(...
          [kpar*ones(numel(P), 1), P],...
          op_pos, mesh_pos, BCstruct_pos,...
          op_neg, mesh_neg, BCstruct_neg,...
          op_int, mesh_int, problem_type,...
          parallel_use, parallel_numcores, plot_coefficients);
    
    bt = bt.refine(indicator_function);

    %% Post-processing
    save(['outputs/EScurves_', num2str(idK)], 'bt');
    % bt.visualize_cache;
  end

  %%
  % figure;
  % for idK = 1:numKpoints
  %   load(['outputs/EScurves_', num2str(idK)]);
  %   bt.visualize_cache;
  %   pause;
  % end

  %%
  figure; hold on;
  for idK = 1:numKpoints
    kpar = kpars(idK);
    load(['outputs/EScurves_', num2str(idK)]);
    [EScurves.XnanBB, EScurves.Xpeaks] = bt.find_nan_boundary_and_peaks(1e5);
    % 
    for idI = 1:size(EScurves.XnanBB, 1)
      plot([kpar, kpar], [EScurves.XnanBB(idI, 1), EScurves.XnanBB(idI, 2)], 'Color', [1 1 1]*100/255);
    end

    for idI = 1:numel(EScurves.Xpeaks)
      plot(kpar, EScurves.Xpeaks(idI), 'r*');
    end
  end

% end

%% Auxiliary functions
function [mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = construct_mesh_and_FEmatrices(pbinputs, HcObj)
  
  %% Coefficients
  % Make sure the edge coefficients are coprime integers
  a1 = HcObj.edge.a1;
  b1 = HcObj.edge.b1;
  [G, x, y] = gcd(a1, b1);  % a1*x + b1*y = G
  a2 = -y; b2 = x;

  if (G ~= 1)
    error('a1 et b1 doivent être premiers entre eux.');
  end

  % Edge vectors
  edge_vec1      =  a1 * HcObj.vecPer1  + b1 * HcObj.vecPer2;
  edge_vec2      =  a2 * HcObj.vecPer1  + b2 * HcObj.vecPer2;
  edge_dual_vec1 =  b2 * HcObj.dualVec1 - a2 * HcObj.dualVec2;
  edge_dual_vec2 = -b1 * HcObj.dualVec1 + a1 * HcObj.dualVec2;

  % Transformation matrices
  Rmat = [edge_vec1, edge_vec2];
  Tmat = [edge_dual_vec1'; edge_dual_vec2'] / (2*pi); % Inverse of Rmat
  Tvec = Tmat' * [1; 0];

  % Potentials
  Vpot = @(x) HcObj.V((Rmat * x(:, 1:2)')');
  Wpot = @(x) HcObj.W((Rmat * x(:, 1:2)')');
  Apot = @(x) HcObj.A((Rmat * x(:, 1:2)')');
  sigma2 = [0 -1i; 1i 0];

  % (+) half-guide
  op_pos.Tmat = Tmat;
  op_pos.vect = Tvec;
  op_pos.funQ = @(x) Vpot(x) + pbinputs.delta * Wpot(x);
  op_pos.funP = @(x) kron(eye(2), ones(size(x, 1), 1)) + pbinputs.delta * kron(sigma2, Apot(x));
  op_pos.funR = 1;

  % (-) half-guide
  op_neg.Tmat = Tmat;
  op_neg.vect = Tvec;
  op_neg.funQ = @(x) Vpot(x) - pbinputs.delta * Wpot(x);
  op_neg.funP = @(x) kron(eye(2), ones(size(x, 1), 1)) - pbinputs.delta * kron(sigma2, Apot(x));
  op_neg.funR = 1;

  % Interior domain
  if strcmpi(pbinputs.problem_type, 'interior')
    yInt = ceil(1/(2*pi*pbinputs.delta));
    BBint = [-yInt, yInt];
    op_int.Tmat = Tmat;
    op_int.vect = Tvec;
    op_int.funQ = @(x) Vpot(x) + pbinputs.delta * HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Wpot(x);
    op_int.funP = @(x) kron(eye(2), ones(size(x, 1), 1)) + pbinputs.delta * kron(sigma2, HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Apot(x));
    op_int.funR = 1;
  else
    op_int = [];
    BBint = [];
  end

  %% Initialization
  [mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = initialize_edge_states_rational(...
    op_pos, pbinputs.numNodesXpos, pbinputs.numNodesYpos, pbinputs.bc_basis_functions_pos,...
    op_neg, pbinputs.numNodesXneg, pbinputs.numNodesYneg, pbinputs.bc_basis_functions_neg,...
    op_int, BBint, pbinputs.numNodesXint, pbinputs.numNodesYint,...
    pbinputs.problem_type...
  );

end
