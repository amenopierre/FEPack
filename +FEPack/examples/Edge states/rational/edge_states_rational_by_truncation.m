%% Compute edge states as eigenvalues of a truncated operator.

%% Import FEPack
clear; clc;
import FEPack.*

%% Problem-related parameters
HcObj = applications.HoneycombObject('none', 'none', 'none');

% Honeycomb lattice potentials
% centers = [-1/sqrt(3), 1/sqrt(3); 0, 0];
% ampsV   = [0;  0];
% ampsW   = [10; 10];
% rads    = [0.2; 0.2];
% HcObj.V = @(x) FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsV, rads);
% HcObj.W = @(x) 1 + FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsW, rads);

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
HcObj.edge.a1 =  1;
HcObj.edge.b1 =  -0.25*sqrt(2);

% Domain wall 
HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);
pbinputs.delta = 0.5;

%% Parallel quasi-momentum range
pbinputs.minKpar = -pi;
pbinputs.maxKpar = +pi;
pbinputs.numKpar = 65;
pbinputs.center_parallel_quasi_momentum_around_K = false;

%% Number of eigenvalues
pbinputs.numEigs = 1e3;

%% Meshes and boundary conditions
% Truncated domain range
pbinputs.domain_range = [-10, 10];

% # nodes per unit for X/Y edge of domain
pbinputs.numNodesX = 64;
pbinputs.numNodesY = 64; 

%% Misc. 
% Should you use parallel computing?
opts.parallel_use = false;
opts.parallel_numcores = 4;

% Plot dispersion functions, save them
opts.plot_disp = false;
opts.save_disp = true;

% [kpars, eigvals] = launch_edge_states_rational_by_truncation(pbinputs, HcObj, opts);

% %% Plot eigenvalues
% figure;
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% hold on;
% for idI = 1:pbinputs.numKpar
%   plot(kpars(idI), real(eigvals(idI, :)), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
% end

%%
% function [kpars, eigvals] = launch_edge_states_rational_by_truncation(pbinputs, HcObj, opts)

  %% Operators
  % Make sure the edge coefficients are coprime integers
  a1 = HcObj.edge.a1;
  b1 = HcObj.edge.b1;
  % [G, x, y] = gcd(a1, b1);  % a1*x + b1*y = G
  % a2 = -y; b2 = x;
  
  % if (G ~= 1)
  %   error('a1 et b1 doivent être premiers entre eux.');
  % end
  a2 = 0; b2 = 1;

  % Edge vectors
  edge_vec1      =  a1 * HcObj.vecPer1  + b1 * HcObj.vecPer2;
  edge_vec2      =  a2 * HcObj.vecPer1  + b2 * HcObj.vecPer2;
  edge_dual_vec1 =  b2 * HcObj.dualVec1 - a2 * HcObj.dualVec2;
  edge_dual_vec2 = -b1 * HcObj.dualVec1 + a1 * HcObj.dualVec2;

  % Transformation matrices
  Rmat = [edge_vec1, edge_vec2];
  Tmat = [edge_dual_vec1'; edge_dual_vec2'] / (2*pi); % Inverse of Rmat
  vect = Tmat' * [1; 0];

  % Potentials
  Vpot = @(x) HcObj.V((Rmat * x(:, 1:2)')');
  Wpot = @(x) HcObj.W((Rmat * x(:, 1:2)')');
  Apot = @(x) HcObj.A((Rmat * x(:, 1:2)')');
  sigma2 = [0 -1i; 1i 0];
  funP = @(x) kron(eye(2), ones(size(x, 1), 1)) + pbinputs.delta * kron(sigma2, HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Apot(x));
  funQ = @(x) Vpot(x) + pbinputs.delta * HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Wpot(x);
  funR = 1;

  %% Mesh
  numNodesY = ceil(abs(pbinputs.domain_range(2) - pbinputs.domain_range(1)) * pbinputs.numNodesY);
  mesh = FEPack.meshes.MeshRectangle(1, [0 1], pbinputs.domain_range, pbinputs.numNodesX, numNodesY);
  dom  = mesh.domain('volumic');
  edge0x = mesh.domain('xmin');
  edge1x = mesh.domain('xmax');
  edgeYmin = mesh.domain('ymin');
  edgeYmax = mesh.domain('ymax');

  %% FE matrices
  u = FEPack.pdes.PDEObject; 
  v = dual(u);

  gradu_gradv = (funP * (Tmat' * grad2(u))) * (Tmat' * grad2(v));
  gradu_vectv = (funP * (Tmat' * grad2(u))) * (vect * v);
  vectu_gradv = (funP * (vect * u)) * (Tmat' * grad2(v));
  vectu_vectv = (funP * (vect * u)) * (vect * v);
  funQ_u_v    = (funQ * u) * v;
  funR_u_v    = (funR * u) * v;

  fprintf('FE matrices computation\n');
  mat_gradu_gradv = FEPack.pdes.Form.intg(dom, gradu_gradv);
  mat_gradu_vectv = FEPack.pdes.Form.intg(dom, gradu_vectv);
  mat_vectu_gradv = FEPack.pdes.Form.intg(dom, vectu_gradv);
  mat_vectu_vectv = FEPack.pdes.Form.intg(dom, vectu_vectv);
  mat_funQ_u_v    = FEPack.pdes.Form.intg(dom, funQ_u_v);
  mat_funR_u_v    = FEPack.pdes.Form.intg(dom, funR_u_v);
  fprintf('FE matrices computation done.\n');

  %% Boundary conditions
  ecs = ((u|edgeYmin) == 0.0) & ((u|edgeYmax) == 0.0)...
      & (((u|edge0x) - (u|edge1x)) == 0.0);
  ecs.applyEcs;
  PP = ecs.P;

  %% Compute eigenvalues
  % kpars = linspace(pbinputs.minKpar, pbinputs.maxKpar, pbinputs.numKpar);
  kpars = HcObj.highSymK' * edge_vec1;
  pbinputs.numKpar = numel(kpars);
  numEigs = pbinputs.numEigs;
  eigvals = zeros(pbinputs.numKpar, numEigs);
  eigfuns = cell(pbinputs.numKpar, 1);

  for idI = 1:pbinputs.numKpar
    fprintf('%d sur %d\n', idI, pbinputs.numKpar);

    kpar   = kpars(idI);

    % FE matrices
    AA = mat_gradu_gradv...
          + 1i * kpar     * mat_vectu_gradv...
          - 1i * kpar     * mat_gradu_vectv...
          + (kpar * kpar) * mat_vectu_vectv...
          + mat_funQ_u_v;

    BB = mat_funR_u_v;

    AA0 = PP * AA * PP';
    BB0 = PP * BB * PP';

    % Solve eigenvalue problem
    [eigfuns{idI}, D] = eigs(AA0, BB0, numEigs, 'smallestabs');
    eigvals(idI, :) = diag(D).';
    
  end


  % % currEigfun = PP' * eigfuns{1}(:, 8);
  % %%
  % perFun = PP' * eigfuns{1}(:, 8);
  % trisurf(mesh.triangles, mesh.points(:, 1), mesh.points(:, 2), real(perFun));
  % axis([0 1 -1 1]); colormap jet;
  % shading interp; view(2); set(gca, 'DataAspectRatio', [1 1 1]);
  % 
  % %%
  % figure;
  % mshCell = FEPack.meshes.MeshRectangle(1, [0 1], [0 1], pbinputs.numNodesX, pbinputs.numNodesY);
  % domCell = mshCell.domain('volumic');
  % 
  % NcellX = 8;
  % NcellY = 1;
  % 
  % VV = cell(2*NcellX, 2*NcellY);
  % 
  % for idX = -NcellX:NcellX-1
  %   for idY = -NcellY:NcellY-1
  %     disp([idX, idY])
  % 
  %     X = mshCell.points(:, 1) + idX;
  %     Y = mshCell.points(:, 2) + idY;
  %     TX = (Tmat * [X, Y].').';
  %     TXmod = [mod(TX(:, 1), 1), TX(:, 2)];
  % 
  %     structLoc = dom.locateInDomain([TXmod, zeros(size(X))]);
  %     elts = dom.elements(structLoc.elements, :)';
  %     elts = elts(:);
  %     coos = structLoc.barycoos';
  %     coos = coos(:);
  % 
  %     exp_k_dot_x = exp(1i * kpars(1) * TX(:, 1));
  %     % V_exp_k_dot_x = exp_k_dot_x .* perFun;
  % 
  %     % VV{idX + NcellX + 1, idY + NcellY + 1} 
  %     val = exp_k_dot_x .* reshape(sum(reshape(coos .* perFun(elts), size(structLoc.barycoos, 2), []), 1), mshCell.numPoints, []);
  %     %
  %     trisurf(mshCell.triangles, X, Y, real(val)); hold on;
  %   end
  % end
  % 
  % shading interp; view(2);
  % set(gca, 'DataAspectRatio', [1 1 1])
  % % trisurf(mesh.triangles, mesh.points(:, 1), mesh.points(:, 2), real())

% end

