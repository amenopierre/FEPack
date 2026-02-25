% function main_edge_states_rational_kpar(a1, b1)
  %% Import FEPack
  import FEPack.*
  % if (nargin < 2), b1 = 0; end
  % if (nargin < 1), a1 = 1; end
  clear; clc;
  a1 = 0; b1 = 1;
  

  %% Problem-related parameters
  HcObj = tools.HoneycombObject('none', 'none', 'none');

  % Honeycomb lattice potentials
  % centers = [-1/sqrt(3), 1/sqrt(3); 0, 0];
  % ampsV   = [0;  0];
  % ampsW   = [10; 10];
  % rads    = [0.2; 0.2];
  % HcObj.V = @(x) FEPack.tools.FUN_atomic_potential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsV, rads);
  % HcObj.W = @(x) 1 + FEPack.tools.FUN_atomic_potential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsW, rads);

  HcObj.V = @(x) 20 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
                 20 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
                 20 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

  HcObj.W = @(x)  1 * sin(x(:, 1:2) *  HcObj.dualVec1) +...
                  1 * sin(x(:, 1:2) *  HcObj.dualVec2) +...
                  1 * sin(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

  HcObj.A = @(x)  0 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
                  0 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
                  0 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

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
  pbinputs.numNodesXpos = 16;
  pbinputs.numNodesYpos = 16; 
  pbinputs.numNodesXneg = 16;
  pbinputs.numNodesYneg = 16;
  pbinputs.numNodesXint = 16;
  pbinputs.numNodesYint = 16;

  % Basis functions for boundary conditions ("Lagrange" or "Fourier")
  pbinputs.bc_basis_functions_pos = "Lagrange";
  pbinputs.bc_basis_functions_neg = "Lagrange";

  %% Misc. 
  % Should you use parallel computing?
  parallel_use = false;
  parallel_numcores = 2;

  % plot the coefficients
  plot_coefficients = false;

  % % Plot dispersion functions, save them
  % opts.plot_disp = false;
  % opts.save_disp = true;

  %% Initialization
  [mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = construct_mesh_and_FEmatrices(pbinputs, HcObj);

  %% Identify bands and gaps
  %  |||||||||||||||||||||||
  numKpoints = 8;
  kpars = linspace(-pi, pi, numKpoints);
  minEgy = pbinputs.minEgy;
  maxEgy = pbinputs.maxEgy;
  problem_type = pbinputs.problem_type;
  
  rootBB = [minEgy, maxEgy];
  maxDepth = 10;
  force_depth = 8;
  force_depth = min(force_depth, maxDepth);

  on_band_border = @(fvals) (any(isnan(fvals), 2) & ~all(isnan(fvals), 2));
  early_split    = @(depth) (depth <= force_depth);

  % The rule: split leaves on the border of continuous spectrum XOR
  % min depth hasn't been reached yet or the indicator function is 
  % above a threshold.
  rules = @(fvals, depth) on_band_border(fvals) | early_split(depth);

  BBband = cell(numKpoints, 1); 
  BBgap = cell(numKpoints, 1);

  figure; hold on;

  for idK = 1:numKpoints
    clc; idK
    kpar = kpars(idK);

    bt = FEPack.meshes.BinaryTree(maxDepth, rootBB, rules);

    indicator_function = @(P) edge_states_rational(...
          [kpar*ones(numel(P), 1), P],...
          op_pos, mesh_pos, BCstruct_pos,...
          op_neg, mesh_neg, BCstruct_neg,...
          op_int, mesh_int, problem_type,...
          parallel_use, parallel_numcores, plot_coefficients, true);
    
    bt = bt.refine(indicator_function);

    % Post-processing
    XnanBB = bt.find_nan_boundary_and_peaks(1e3);

    if (size(XnanBB, 1) == 0)
      BBgap{idK} = [minEgy, maxEgy];
    else
      BBgap{idK}(1, :) = [minEgy, XnanBB(1, 1)];

      for idB = 2:size(XnanBB, 1)
        BBgap{idK}(idB, :) = [XnanBB(idB-1, 2) XnanBB(idB, 1)];
      end

      % for idI = 1:size(XnanBB, 1)
      %   plot([kpar, kpar], [XnanBB(idI, 1), XnanBB(idI, 2)], 'Color', [1 1 1]*100/255);
      % end
    end
    BBband{idK} = XnanBB;
  end 
  
  %% Compute eigenvalues within gaps
  %  |||||||||||||||||||||||||||||||
  tol = 1e-6;
  
  Eval = cell(idK, 1);

  for idK = 1:numKpoints
    clc; idK
    kpar = kpars(idK);
    gapsK = BBgap{idK};
    Eval{idK} = cell(size(gapsK, 1), 1);

    for idB = 1:size(gapsK, 1)

      Ngap = ceil(64 * abs(gapsK(idB, 2) - gapsK(idB, 1))) + 2;
      Epts = linspace(gapsK(idB, 1), gapsK(idB, 2), Ngap).';
      Epts(1:2) = [];
      Epts(end-1:end) = [];

      [~, eigvals] = edge_states_rational(...
          [kpar*ones(numel(Epts), 1), Epts],...
          op_pos, mesh_pos, BCstruct_pos,...
          op_neg, mesh_neg, BCstruct_neg,...
          op_int, mesh_int, problem_type,...
          parallel_use, parallel_numcores, plot_coefficients);

      eigvals = cell2mat(eigvals);

      esId = find( abs(eigvals - Epts*ones(1, size(eigvals, 2))) < tol);
      [esId, ~] = ind2sub(size(eigvals), esId);
      Eval{idK}{idB} = Epts(esId);
    end
  end

  %
  for idK = 1:numKpoints
    kpar = kpars(idK);
    gapsK = BBgap{idK};

    for idB = 1:size(gapsK, 1)
      plot(kpar*ones(numel(Eval{idK}{idB}), 1), Eval{idK}{idB}, 'b*'); hold on;
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
    op_int.funP = @(x) kron(eye(2), ones(size(x, 1), 1)) - pbinputs.delta * kron(sigma2, HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Apot(x));
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
