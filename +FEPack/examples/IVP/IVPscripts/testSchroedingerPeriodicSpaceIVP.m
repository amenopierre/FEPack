%% Import FEPack
import FEPack.*
clear; clc;
close all;

%  ***************************** %
%% Generate mesh and FE matrices %
%  ***************************** %
honeycomb_mesh_FEmatrices;

folder_name = '/home/pierre/Documents/Recherche/Code/FEPack/+FEPack/outputs/';
system(['mkdir ', folder_name, 'IVPout']);
system(['rm ',  folder_name, 'IVPout/*']); % In case the folder already exists

%  ***************************************** %
%% Time-dependent variables and initial data %
%  ***************************************** %
theta = 0.5;
delta_t = 1e-3;
final_t = 1e-1;

% Initial data
initial_data_before_transformation = @(x) FEPack.tools.cutoff(sqrt((x(:, 1)-0.5).^2 + x(:, 2).^2), -0.5, 0.5);
initial_data = @(x) initial_data_before_transformation((pbinputs.Rmat * x(:, 1:2)')');

% asol0 = 30;
% initial_data = @(x) exp(-asol0 * x(:, 2).^2);% cutoff(x, -0.5, 0.5);

% vecInitial = initial_data(mshInt.points);
rhsInt = @(x, t) zeros(size(x, 1), 1);

% Floquet dual variable
numFBpoints = 64;
FBpoints = linspace(-pi, pi, numFBpoints);

%  ************ %
%% Floquet loop %
%  ************ %
parfor idFB = 1:numFBpoints

  % Current Floquet variable
  kpar = FBpoints(idFB);

  %% Floquet-Bloch transform of initial data
  vecInitial = FEPack.tools.BlochTransform(mshInt.points, kpar, initial_data, 1);

  %% Floquet-Bloch transform of initial data
  rhsIntFB = @(t) FEPack.tools.BlochTransform(mshInt.points, kpar, (@(x) rhsInt(x, t)), 1);

  rhsTheta = @(idT) op_int.mat_funR_u_v * (...
              theta  * rhsIntFB((idT  ) * delta_t) +...
        (-1 + theta) * rhsIntFB((idT-1) * delta_t));

  %% Construct FE matrices
  HHpos = op_pos.mat_gradu_gradv...
    + 1i * kpar     * op_pos.mat_vectu_gradv...
    - 1i * kpar     * op_pos.mat_gradu_vectv...
    + (kpar * kpar) * op_pos.mat_vectu_vectv...
    + op_pos.mat_funQ_u_v;

  AApos =       theta  * HHpos - (1i/delta_t) * op_pos.mat_funR_u_v;
  BBpos = (-1 + theta) * HHpos - (1i/delta_t) * op_pos.mat_funR_u_v;

  HHneg = op_neg.mat_gradu_gradv...
    + 1i * kpar     * op_neg.mat_vectu_gradv...
    - 1i * kpar     * op_neg.mat_gradu_vectv...
    + (kpar * kpar) * op_neg.mat_vectu_vectv...
    + op_neg.mat_funQ_u_v;

  AAneg =       theta  * HHneg - (1i/delta_t) * op_neg.mat_funR_u_v;
  BBneg = (-1 + theta) * HHneg - (1i/delta_t) * op_neg.mat_funR_u_v;

  HHint = op_int.mat_gradu_gradv...
    + 1i * kpar     * op_int.mat_vectu_gradv...
    - 1i * kpar     * op_int.mat_gradu_vectv...
    + (kpar * kpar) * op_int.mat_vectu_vectv...
    + op_int.mat_funQ_u_v;

  AAint =       theta  * HHint - (1i/delta_t) * op_int.mat_funR_u_v;
  BBint = (-1 + theta) * HHint - (1i/delta_t) * op_int.mat_funR_u_v;

  outFB = schroedingerPeriodicGuideIVP(...
    delta_t, final_t, 2,...
    AApos, BBpos, mshPos, op_pos.BCstruct, pbinputs.numCellsPos,...
    AAneg, BBneg, mshNeg, op_neg.BCstruct, pbinputs.numCellsNeg,...
    AAint, BBint, vecInitial, mshInt, op_int.spBint_pos, op_int.spBint_neg,...
    rhsTheta);

  FEPack.tools.parsave([folder_name, 'IVPout/outFB_', num2str(idFB)], outFB, true);

end

% clear outFB;

%  ************************* %
%% Inverse Floquet transform %
%  ************************* %
numCellsInfinite = 3;
numTsteps = ceil(final_t / delta_t);

%% Positive side
op_pos.vecCells = [2*numCellsInfinite pbinputs.numCellsPos]; NuPos = prod(op_pos.vecCells);
op_neg.vecCells = [2*numCellsInfinite pbinputs.numCellsNeg]; NuNeg = prod(op_neg.vecCells);
op_int.vecCells = [2*numCellsInfinite 1]; NuInt = prod(op_int.vecCells);

[I1, I2] = ind2sub(op_pos.vecCells, 1:NuPos); op_pos.pointsIds = [I1; I2]; % 2-by-NuPos
[I1, I2] = ind2sub(op_neg.vecCells, 1:NuNeg); op_neg.pointsIds = [I1; I2]; % 2-by-NuNeg
[I1, I2] = ind2sub(op_int.vecCells, 1:NuInt); op_int.pointsIds = [I1; I2]; % 2-by-NuInt

tauPos = op_pos.pointsIds(1, :) - numCellsInfinite' * ones(1, NuPos) - 1; % 1-by-NuPos
tauNeg = op_neg.pointsIds(1, :) - numCellsInfinite' * ones(1, NuNeg) - 1; % 1-by-NuNeg
tauInt = op_int.pointsIds(1, :) - numCellsInfinite' * ones(1, NuInt) - 1; % 1-by-NuInt

W = 2*pi ./ (numFBpoints - 1); % 1-by-1
W = W * sqrt(1 / (2*pi));
out.pos = zeros(mshPos.numPoints, NuPos, numTsteps);
out.neg = zeros(mshNeg.numPoints, NuNeg, numTsteps);
out.int = zeros(mshInt.numPoints, NuInt, numTsteps);


% val = zeros(3, 2, 8);
% parfor i = 1:8
%   val(:, :, i) = val(:, :, i) + 5;
% end

for idFB = 1:numFBpoints
  % Load Floquet-Bloch transform
  kpar = FBpoints(idFB);
  outFB = load([folder_name, 'IVPout/outFB_', num2str(idFB)]);

  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangle rule.

  %% Positive side
  exp_k_dot_x = exp(1i * mshPos.points(:, 1) * kpar); % Npos-by-1
  exp_k_dot_tau = exp(1i * kpar * tauPos); % 1-by-NuPos
  
  for idT = 1:numTsteps
    U_TFB = outFB.pos.sol(:, op_pos.pointsIds(2, :), idT); % Npos-by-NuPos
    out.pos(:, :, idT) = out.pos(:, :, idT) + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
  end

  %% Negative side
  exp_k_dot_x = exp(1i * mshNeg.points(:, 1) * kpar); % Nneg-by-1
  exp_k_dot_tau = exp(1i * kpar * tauNeg); % 1-by-NuNeg
  
  for idT = 1:numTsteps
    U_TFB = outFB.neg.sol(:, op_neg.pointsIds(2, :), idT); % Nneg-by-NuNeg
    out.neg(:, :, idT) = out.neg(:, :, idT) + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
  end

  %% Interior domain
  exp_k_dot_x = exp(1i * mshInt.points(:, 1) * kpar); % Nint-by-1
  exp_k_dot_tau = exp(1i * kpar * tauInt); % 1-by-NuInt
  
  for idT = 1:numTsteps
    U_TFB = outFB.int.sol(:, idT) * ones(1, NuInt); % Nint-by-NuInt
    out.int(:, :, idT) = out.int(:, :, idT) + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
  end
end

% error;
%% Plot solution
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idT = 1:numTsteps
  for idI = 1:op_pos.vecCells(1)
    % Positive side
    for idS = 1:op_pos.vecCells(2)
      idcell = sub2ind(op_pos.vecCells, idI, idS);

      X = mshPos.points(:, 1) + (idI - numCellsInfinite - 1);
      Y = mshPos.points(:, 2) + (idS - 1) + pbinputs.IDb(2);

      trisurf(mshPos.triangles, X, Y, real(out.pos(:, idcell, idT)));
      hold on;
    end
    % shading interp; view(2); 
    % pause;

    % Negative side
    for idS = 1:op_neg.vecCells(2)
      idcell = sub2ind(op_neg.vecCells, idI, idS);

      X = mshNeg.points(:, 1) + (idI - numCellsInfinite - 1);
      Y = mshNeg.points(:, 2) - (idS - 1) + pbinputs.IDb(1);

      trisurf(mshNeg.triangles, X, Y, real(out.neg(:, idcell, idT)));
      hold on;
    end
    % shading interp; view(2); 
    % pause;

    % Interior domain
    X = mshInt.points(:, 1) + (idI - numCellsInfinite - 1);
    Y = mshInt.points(:, 2);

    trisurf(mshInt.triangles, X, Y, real(out.int(:, idI, idT)));
    hold on;
    % shading interp; view(2); 
    % pause;
  end

  % view(90, 2);
  view(2);
  shading interp; 
  colorbar('TickLabelInterpreter', 'latex');
  colormap jet;
  set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);

  if (idT == 1)
    pause;
  else
    pause(0.05);
  end
end


