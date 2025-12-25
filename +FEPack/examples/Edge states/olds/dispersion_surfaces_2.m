% Compute dispersion curves associated
% to a honeycomb lattice potential
clear; clc;
import FEPack.*

% Lattice vectors
v1 = [0.5*sqrt(3);  0.5];
v2 = [0.5*sqrt(3); -0.5];

% Dual lattice vectors
k1 = (4*pi/sqrt(3)) * [0.5;  0.5*sqrt(3)];
k2 = (4*pi/sqrt(3)) * [0.5; -0.5*sqrt(3)];
Khs = (k1 - k2) / 3; % High-symmetry quasi-momentum
brillouinVerts = [Khs, -Khs + k1, Khs + k2, -Khs, Khs - k1, -Khs - k2];

% Honeycomb lattice potential
A = [0; 0];
B = [1/sqrt(3); 0];
ampA = 20;
ampB = 0;   % Attention : V n'est pas paire si ampB ≠ 0
% Vpot = @(x) zeros(size(x, 1), 1);
Vpot = @(x) atomicPotential(x, v1, v2, A, B, ampA, ampB);
% Vpot = @(x) opticalPotential(x, k1, k2, true);
% Vpot = @(x) 10*trigonometricPolynomial(x, k1, k2, true);

Wpot = @(x) opticalPotential(x, k1, k2, false);
% Wpot = @(x) trigonometricPolynomial(x, k1, k2, false);
delta = 0;
Qpot = @(x) Vpot(x) + delta * Wpot(x);

numNodesX = 100;
meshXY = meshes.MeshRectangle(1, [0 1], [0 1], numNodesX, numNodesX);
cellXY = meshXY.domain('volumic');
edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

% %% Represent coefficient
% numX = 2;
% numY = 2;
% R = [v1, v2];
% transQpot = @(x) Qpot((R * x(:, 1:2)')');
% figure;
% for idX = -numX:numX-1
%   for idY = -numY:numY-1
%     X = meshXY.points(:, 1) + idX;
%     Y = meshXY.points(:, 2) + idY;

%     % trisurf(meshXY.triangles, X, Y, Qpot([X, Y]));
%     trisurf(meshXY.triangles, X, Y, transQpot([X, Y])); 
%     hold on;
%     shading interp;
%     view(2);
%     set(gca, 'DataAspectRatio', [1 1 1]);
%   end
% end
% error();
% %%

% Number of eigenvalues to compute
numEVP = 4;

% Forms
u = pdes.PDEObject; v = dual(u);
R = [v1, v2];
T = R \ eye(2);
Te1 = T' * [1; 0];
Te2 = T' * [0; 1];
transQpot = @(x) Qpot((R * x(:, 1:2)')');

fprintf('Calcul des matrices EF\n');
mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, (T' * grad2(u)) * (T' * grad2(v)));
mat_gradu_vec1v = FEPack.pdes.Form.intg(cellXY, (T' * grad2(u)) * (Te1 * v));
mat_gradu_vec2v = FEPack.pdes.Form.intg(cellXY, (T' * grad2(u)) * (Te2 * v));
mat_vec1u_gradv = FEPack.pdes.Form.intg(cellXY, (Te1 * u) * (T' * grad2(v)));
mat_vec2u_gradv = FEPack.pdes.Form.intg(cellXY, (Te2 * u) * (T' * grad2(v)));
mat_vec1u_vec1v = FEPack.pdes.Form.intg(cellXY, (Te1 * u) * (Te1 * v));
mat_vec1u_vec2v = FEPack.pdes.Form.intg(cellXY, (Te1 * u) * (Te2 * v));
mat_vec2u_vec1v = FEPack.pdes.Form.intg(cellXY, (Te2 * u) * (Te1 * v));
mat_vec2u_vec2v = FEPack.pdes.Form.intg(cellXY, (Te2 * u) * (Te2 * v));
mat_Qpot_u_v    = FEPack.pdes.Form.intg(cellXY, (transQpot * u) * v);
mat_u_v = FEPack.pdes.Form.intg(cellXY, u * v);
fprintf('Calcul terminé.\n');

% Periodic boundary conditions
ecs = (((u|edge0x) - (u|edge1x)) == 0.0) & (((u|edge0y) - (u|edge1y)) == 0.0);
ecs.applyEcs;

%% Compute the band slice
numNodesK = 6;
numNodesK1 = 32;% ceil((4*pi/3) * numNodesK);
numNodesK2 = 32;% ceil((2*pi/sqrt(3)) * numNodesK);
% meshK = meshes.MeshRectangle(1, [-k1(1) k1(1)], [-Khs(2) Khs(2)], numNodesK1, numNodesK2);

numSubK = 3;
subKx = linspace(-(k1(1)+0.5),   (k1(1)+0.5), numSubK);
subKy = linspace(-(Khs(2)+0.5), (Khs(2)+0.5), numSubK);

for idPx = 1:numSubK-1
  for idPy = 1:numSubK-1

    meshK = meshes.MeshRectangle(1, [subKx(idPx) subKx(idPx+1)], [subKy(idPy) subKy(idPy+1)], numNodesK1, numNodesK2);
    meshK.MeshRectangleLocallyRefined(0, [subKx(idPx) subKx(idPx+1)], [subKy(idPy) subKy(idPy+1)], numNodesK1, numNodesK2, brillouinVerts, 0.01*ones(1, 6));

    Kpts = meshK.points(:, 1:2);
    numK = size(Kpts, 1);

    eigvals = zeros(numK, numEVP);

    P = ecs.P;

    for idI = 1:numK
      fprintf('%d sur %d\n', idI, numK);
      
      kv = Kpts(idI, :);
      kx = kv * v1;
      ky = kv * v2;

      % FE matrix
      AA = mat_gradu_gradv +...
        + 1i * kx * mat_vec1u_gradv...
        + 1i * ky * mat_vec2u_gradv...
        - 1i * kx * mat_gradu_vec1v...
        - 1i * ky * mat_gradu_vec2v...
        + (kx * kx) * mat_vec1u_vec1v...
        + (kx * ky) * mat_vec1u_vec2v...
        + (ky * kx) * mat_vec2u_vec1v...
        + (ky * ky) * mat_vec2u_vec2v...
        + mat_Qpot_u_v;

      BB = mat_u_v;

      AA0 = P * AA * P';
      BB0 = P * BB * P';

      % Solve eigenvalue problem
      [~, D] = eigs(AA0, BB0, numEVP, 'smallestabs');
      eigvals(idI, :) = diag(D).';
    end

    eigvals = real(eigvals);

    %%
    % figure;
    % % hold off;
    % % plot(lambdas, 0, 'rx');
    % for idJ = 1:numEVP

    %   trisurf(meshK.triangles, Kpts(:, 1), Kpts(:, 2), eigvals(:, idJ)); hold on;

    % end
    % shading interp

    % Edirac = 18.35;
    % zlim([Edirac-0.05, Edirac+0.05]);
    % clim([Edirac-0.05, Edirac+0.05]);
    set(gca, 'DataAspectRatio', [1 1 1]);
    % view(2);
    % ylim([0 30])

    %%
    fid = fopen(['outputs/D0_Kdomaine_X', int2str(idPx), '_Y', int2str(idPy), '.txt'], 'w+');
    % Kdom = [];
    % for idI = 1:2:meshK.numTriangles;
    %   rect = unique(meshK.triangles(idI:idI+1, :));
    %   Kdom = [Kdom; rect(:)'-1];
    % end
    fprintf(fid, 'pt1\t pt2\t pt3\n');
    fprintf(fid, '%0.6d\t %0.6d\t %0.6d\n', -1+meshK.triangles.');
    fclose(fid);

    fid = fopen(['outputs/D0_dispersion_surfaces_X', int2str(idPx), '_Y', int2str(idPy), '.txt'], 'w+');
    fprintf(fid, 'x\t y');
    for idJ = 1:numEVP
      fprintf(fid, ['\t E', int2str(idJ)]);
    end
    fprintf(fid, '\n');
    fprintf(fid, ['%0.3e\t%0.3e', repmat('\t%0.3e', 1, numEVP), '\n'], [Kpts, eigvals].');
    fclose(fid);

  end
end

kv = Khs.';
kx = kv * v1;
ky = kv * v2;

% FE matrix
AA = mat_gradu_gradv +...
  + 1i * kx * mat_vec1u_gradv...
  + 1i * ky * mat_vec2u_gradv...
  - 1i * kx * mat_gradu_vec1v...
  - 1i * ky * mat_gradu_vec2v...
  + (kx * kx) * mat_vec1u_vec1v...
  + (kx * ky) * mat_vec1u_vec2v...
  + (ky * kx) * mat_vec2u_vec1v...
  + (ky * ky) * mat_vec2u_vec2v...
  + mat_Qpot_u_v;

BB = mat_u_v;

AA0 = P * AA * P';
BB0 = P * BB * P';

% Solve eigenvalue problem
[~, D] = eigs(AA0, BB0, numEVP, 'smallestabs');

function val = atomicPotential(x, v1, v2, A, B, ampA, ampB)

  v1 = v1(:);
  v2 = v2(:);
  A = A(:);
  B = B(:);

  R = [v1, v2];
  T = R \ eye(2);
  normR = norm(R, 'fro');

  XmodA = (R * (mod(T * (x(:, 1:2).' - A) + 0.5, 1) - 0.5)).';
  XmodB = (R * (mod(T * (x(:, 1:2).' - B) + 0.5, 1) - 0.5)).';

  val = ampA * FEPack.tools.cutoff(sqrt(XmodA(:, 1).^2 + XmodA(:, 2).^2), -0.4*normR, 0.4*normR) +...
        ampB * FEPack.tools.cutoff(sqrt(XmodB(:, 1).^2 + XmodB(:, 2).^2), -0.2*normR, 0.2*normR);

end

function val = opticalPotential(x, k1, k2, is_even)
  
  if (is_even)
    Xmod = cos(x(:, 1:2) * k1) + cos(x(:, 1:2) * k2) + cos(x(:, 1:2) * (k1 + k2));

    val = FEPack.tools.cutoff(Xmod, -1.5, 1.5);
  else
    Xmod = sin(x(:, 1:2) * k1) + sin(x(:, 1:2) * k2) + sin(x(:, 1:2) * (k1 + k2));

    val = tanh(Xmod);
  end

end

function val = trigonometricPolynomial(x, k1, k2, is_even)
  
  if (is_even)
    val = cos(x(:, 1:2) * k1) + cos(x(:, 1:2) * k2) + cos(x(:, 1:2) * (k1 + k2));
  else
    val = sin(x(:, 1:2) * k1) + sin(x(:, 1:2) * k2) + sin(x(:, 1:2) * (k1 + k2));
  end

end

function val = modBrillouinZone(k)

  % k is 2 by N matrix

  v1 = [0.5*sqrt(3);  0.5];
  v2 = [0.5*sqrt(3); -0.5];

  k1 = (4*pi/sqrt(3)) * [0.5;  0.5*sqrt(3)];
  k2 = (4*pi/sqrt(3)) * [0.5; -0.5*sqrt(3)];

  Khs = (k1 - k2) / 3;

  % Project k onto the unit parallelogram spanned by k1 and k2
  kpar = k1 * mod(v1' * k / (2*pi), 1) + k2 * mod(v2' * k / (2*pi), 1);

  % Find which part of the parallelogram k is in
  pol1 = [[0; 0], 0.5*k2, Khs + k2, -Khs + k1, 0.5*k1];
  pol2 = [k1 + k2, k1 + k2 - 0.5*k2, -Khs + k1, Khs + k2, k1 + k2 - 0.5*k1];
  pol3 = [Khs + k2, 0.5*k2, k2, k1 + k2 - 0.5*k1];
  pol4 = [k1, 0.5*k1, -Khs + k1, k1 + k2 - 0.5*k2];

  % Plot the polygons
  % plot([0, k1(1)], [0, k1(2)], 'k'); hold on;
  % plot([0, k2(1)], [0, k2(2)], 'k');
  % plot([k1(1), k1(1) + k2(1)], [k1(2), k1(2) + k2(2)], 'k');
  % plot([k2(1), k1(1) + k2(1)], [k2(2), k1(2) + k2(2)], 'k');
  % plot(pol1(1, :), pol1(2, :), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
  % plot(pol2(1, :), pol2(2, :), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
  % plot(pol3(1, :), pol3(2, :), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
  % plot(pol4(1, :), pol4(2, :), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
  % set(gca, 'DataAspectRatio', [1 1 1])

  I1 = inpolygon(kpar(1, :), kpar(2, :), pol1(1, :), pol1(2, :));
  I2 = inpolygon(kpar(1, :), kpar(2, :), pol2(1, :), pol2(2, :));
  I3 = inpolygon(kpar(1, :), kpar(2, :), pol3(1, :), pol3(2, :));
  I4 = inpolygon(kpar(1, :), kpar(2, :), pol4(1, :), pol4(2, :));

  val = zeros(size(k));
  val(:, I1) = kpar(:, I1);
  val(:, I2) = kpar(:, I2) - (k1 + k2);
  val(:, I3) = kpar(:, I3) - k2;
  val(:, I4) = kpar(:, I4) - k1;

end
