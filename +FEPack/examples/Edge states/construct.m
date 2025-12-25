clear; clc
out = load('outputs/qt_1_0.mat');
qt = out.qt;
%
H = figure;
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% qt.visualize_cache(@(x) log(abs(x)));
qt.visualize_cache;

% axis([-pi, pi, 7, 16]);
% colorbar;

set(H, 'Position', get(0, 'Screensize'), 'visible', 'off');
% clim([1, 2]);
axis off;
colorbar off;
view(2);
print(H, 'outputs/qt_1_0', '-dpng');
system(['convert outputs/qt_1_0.png -trim outputs/qt_1_0.png']);
close(H);


% subplot(1, 2, 2);
% qt.visualize;%([], 8)% (rootBB);

%%
numKpoints = 256;
kpars = linspace(-pi, pi, numKpoints);
hold on;

for idK = 1:numKpoints
  kpar = kpars(idK);

  EScurves = load(['outputs/out/outputs/EScurves_', num2str(idK)]);

  for idI = 1:size(EScurves.XnanBB, 1)
    plot([kpar, kpar], [EScurves.XnanBB(idI, 1), EScurves.XnanBB(idI, 2)], 'k-x');
  end

  for idI = 1:numel(EScurves.Xpeaks)
    plot(kpar, EScurves.Xpeaks(idI), 'r*');
  end
end

ylim(11.7341 + [-10 10]);





%%
numKpoints = 256;
kpars = linspace(-pi, pi, numKpoints);
figure;
hold on;

fid1 = fopen('outputs/bulk_spectrum.txt', 'w+');
fid2 = fopen('outputs/edge_spectrum.txt', 'w+');
N = 2;
msg = ['kpar\tminB1\tmaxB1' strjoin(arrayfun(@(x) sprintf('\tminB%d\tmaxB%d', x, x), 2:N, 'UniformOutput', false), '') '\n'];
fprintf(fid1, msg);
fprintf(fid2, 'kpar\tEigval');

for idK = 1:numKpoints
  kpar = kpars(idK);

  load(['outputs/out/outputs/EScurves_', num2str(idK)]);
  
  % bt.visualize_cache;
  % pause;
  
  [EScurves.XnanBB, EScurves.Xpeaks] = bt.find_nan_boundary_and_peaks(5e6, 1e2);
  % pause;
  % 

  Ycol = EScurves.XnanBB.';
  Ycol = Ycol(:).';
  fprintf(fid1, ['%0.5e', repmat('\t%0.5e', 1, numel(Ycol)), '\n'], [kpar, Ycol].');

  % disp(size(EScurves.Xpeaks));
  
  if ~isempty(EScurves.Xpeaks)
    fprintf(fid2, '\n%0.5e\t%0.5e', [kpar, EScurves.Xpeaks(1)]);
    
    
    for idI = 2:numel(EScurves.Xpeaks)
      % disp([kpar, EScurves.Xpeaks(idI)])
      fprintf(fid2, '\n%0.5e\t%0.5e', [kpar, EScurves.Xpeaks(idI)]);
    end
  end

  for idI = 1:size(EScurves.XnanBB, 1)
    plot([kpar, kpar], [EScurves.XnanBB(idI, 1), EScurves.XnanBB(idI, 2)], 'Color', [1 1 1]*100/255);
  end

  for idI = 1:numel(EScurves.Xpeaks)
    plot(kpar, EScurves.Xpeaks(idI), 'ro');
  end
end

fclose(fid1);
fclose(fid2);

% axis([-pi, pi, -8, 12])
% axis([-pi, pi, 11.7341 + [-10 10]]);
% set(gca, 'DataAspectRatio', [1 1 1])















%%
numKcell = 6;
numEcell = 4;

folder_name = '/';
% folder_name = 'out_1_0/';
% folder_name = 'out_1_1/';
% folder_name = 'out_2_1/';
% folder_name = 'out_3_1/';
% folder_name = 'out_4_1/';

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

for idK = 1:numKcell
  for idE = 1:numEcell
    try
      out = load(['outputs/' folder_name, 'disp_', num2str(idK), '_', num2str(idE)]);

      % val = abs(out.val);
      val = log(abs(out.val));
      val(isinf(val)) = -1;
      trisurf(out.dispmesh.triangles, out.dispmesh.points(:, 1), out.dispmesh.points(:, 2), val); hold on;

      % trisurf(out.dispmesh.triangles, out.dispmesh.points(:, 1), out.dispmesh.points(:, 2), log(abs(out.val))); hold on;
    catch ME
    end
  end
end

% axis([kpars(1) kpars(end) Egies(1) Egies(end)]);
shading interp;
colormap parula;
view(2);
set(gca, 'FontSize', 16);
colorbar('TickLabelInterpreter', 'latex');
xlabel('$k_\parallel$');
ylabel('$E$');
