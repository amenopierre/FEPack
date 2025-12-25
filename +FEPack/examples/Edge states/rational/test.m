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
ampA = 1;
ampB = 1;   % Attention : V n'est pas paire si ampB ≠ 0
% Vpot = @(x) zeros(size(x, 1), 1);
Vpot = @(x) atomicPotential(x, v1, v2, A, B, ampA, ampB);
% Vpot = @(x) opticalPotential(x, k1, k2, true);
% Vpot = @(x) 10*trigonometricPolynomial(x, k1, k2, true);

% Vpot = @(x) 1 * cos(x(:, 1:2) *  k1) +...
%             1 * cos(x(:, 1:2) *  k2) +...
%             1 * cos(x(:, 1:2) * (k1  + k2));

Wpot = @(x) atomicPotential(x, v1, v2, A, B, 1, -1);
% Wpot = @(x) opticalPotential(x, k1, k2, false);
% Wpot = @(x) trigonometricPolynomial(x, k1, k2, false);
% Wpot = @(x) 1 * sin(x(:, 1:2) *  k1) +...
%             1 * sin(x(:, 1:2) *  k2) +...
%             1 * sin(x(:, 1:2) * (k1  + k2));
delta = 0.5;

edgeslope = -0.5*sqrt(2);
k2edge = -edgeslope * k1 + k2;
kappa = @(y) tanh(1e2*y);

Qpot = @(x) Vpot(x) + delta * kappa(delta * x(:, 1:2) * k2edge) .* Wpot(x);

numNodesX = 100;
meshXY = meshes.MeshRectangle(1, [0 1], [0 1], numNodesX, numNodesX);
cellXY = meshXY.domain('volumic');
edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

%% Represent coefficient
numX = 8;
numY = 8;
R = [v1, v2];
transQpot = @(x) Qpot((R * x(:, 1:2)')');
H = figure;
for idX = -numX:numX-1
  for idY = -numY:numY-1
    X = meshXY.points(:, 1) + idX;
    Y = meshXY.points(:, 2) + idY;

    trisurf(meshXY.triangles, X, Y, Qpot([X, Y]));
    % trisurf(meshXY.triangles, X, Y, transQpot([X, Y])); 
    hold on;
    shading interp;
    view(2);
    set(gca, 'DataAspectRatio', [1 1 1]);
  end
end

set(H, 'Position', get(0, 'Screensize'), 'visible', 'off');
% clim([1, 2]);
axis off;
colorbar off;
view(2);
print(H, 'outputs/medium', '-dpng');
system(['convert outputs/medium.png -trim outputs/medium.png']);
close(H);

error();

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

  val = ampA * FEPack.tools.cutoff(sqrt(XmodA(:, 1).^2 + XmodA(:, 2).^2), -0.3*normR, 0.3*normR) +...
        ampB * FEPack.tools.cutoff(sqrt(XmodB(:, 1).^2 + XmodB(:, 2).^2), -0.3*normR, 0.3*normR);

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
