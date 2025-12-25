%% Import FEPack
import FEPack.*
clear; clc;
close all;

%% Generate mesh and FE matrices
honeycomb_mesh_FEmatrices;

%% Time-dependent variables and initial data
theta = 0.5;
delta_t = 1e-3;
final_t = 1e-1;

% Initial data
initial_data_before_transformation = @(x) FEPack.tools.cutoff(sqrt((x(:, 1)-0.5).^2 + x(:, 2).^2), -0.5, 0.5);
initial_data = @(x) initial_data_before_transformation((Rmat * x(:, 1:2)')');

% asol0 = 30;
% initial_data = @(x) exp(-asol0 * x(:, 2).^2);% cutoff(x, -0.5, 0.5);

op_int.vecInitial = initial_data(mshInt.points);
op_int.rhsInt = @(x, t) zeros(size(x, 1), 1);

% kpar
kpar = 0;

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

%% Right-hand side
rhsTheta = @(idT) op_int.mat_funR_u_v * (...
        theta  * op_int.rhsInt(mshInt.points, (idT  ) * delta_t) +...
  (-1 + theta) * op_int.rhsInt(mshInt.points, (idT-1) * delta_t));

out = schroedingerPeriodicGuideIVP(...
  delta_t, final_t, 2,...
  AApos, BBpos, mshPos, op_pos.BCstruct, pbinputs.numCellsPos,...
  AAneg, BBneg, mshNeg, op_neg.BCstruct, pbinputs.numCellsNeg,...
  AAint, BBint, op_int.vecInitial, mshInt, op_int.spBint_pos, op_int.spBint_neg,...
  rhsTheta);

% error;
%% Figure
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

numTsteps = ceil(final_t / delta_t);

for idT = 1:numTsteps
  for idI = 1:pbinputs.numCellsPos
    X = mshPos.points(:, 1);
    Y = mshPos.points(:, 2) + pbinputs.IDb(2) + idI-1;

    trisurf(mshPos.triangles, X, Y, real(out.pos.sol(:, idI, idT)));
    % trisurf(mshPos.triangles, X, Y, real(initial_data([X, Y])));
    hold on;
  end

  for idI = 1:pbinputs.numCellsNeg
    X = mshNeg.points(:, 1);
    Y = mshNeg.points(:, 2) + pbinputs.IDb(1) - (idI-1);

    trisurf(mshNeg.triangles, X, Y, real(out.neg.sol(:, idI, idT)));
    % trisurf(mshNeg.triangles, X, Y, real(initial_data([X, Y])));
    hold on;
  end

  X = mshInt.points(:, 1);
  Y = mshInt.points(:, 2);
  trisurf(mshInt.triangles, X, Y, real(out.int.sol(:, idT)));
  % trisurf(mshInt.triangles, X, Y, real(initial_data([X, Y])));
  hold off;
  shading interp;
  view(2);
  % view(90, 2);
  set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
  colorbar('TickLabelInterpreter', 'latex');
  colormap jet;
  % clim([-1e-1, 1]);

  if (idT == 1)
    pause;
  else
    pause(0.05);
  end
  % pause;
end


