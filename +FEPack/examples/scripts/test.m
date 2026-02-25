clear; clc;
%%
import FEPack.*

% profile ON;
numNodes2D = 50;
m = load(['pregenMeshes/2D/unstruct_mesh_2D_', int2str(numNodes2D), '_positive.mat']);
mesh2Dpos = m.mesh;
vecperA = [1; 0];
vecperB = [-0.5*sqrt(2); 1];

numCellsSemiInfinite_pos = 3;
numCellsInfinite = 2;

tic;
% mu2D = generate_coefficient('cutoff2', [1; 0], [sqrt(2), 1]);
mu2D = @(x) tools.FUN_per_cutoff_circle(x, vecperA, vecperB, [0.25, 0.75], [-0.1, 0.1]) +...
            tools.FUN_per_cutoff_cuboid(x, vecperA, vecperB, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
% mu2D = @(x) tools.FUN_per_cutoff_cuboid(x, vecperA, vecperB, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
toc;

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


for idS = 1:(2*numCellsSemiInfinite_pos)
  for idI = 1:(2*numCellsInfinite)
    X = mesh2Dpos.points(:, 1) + (idS - numCellsSemiInfinite_pos - 1);
    Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);
    
    trisurf(mesh2Dpos.triangles, X, Y, mu2D([X, Y]));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

for idK = -5:5
  for idL = -5:5
    P = idK * vecperA + idL * vecperB;
    
    plot([P(1), P(1) + vecperA(1)], [P(2), P(2) + vecperA(2)], 'w');
    plot([P(1), P(1) + vecperB(1)], [P(2), P(2) + vecperB(2)], 'w');

  end
end

% profile viewer
xlim([-numCellsSemiInfinite_pos, numCellsSemiInfinite_pos]);
ylim([-numCellsInfinite, numCellsInfinite]);

