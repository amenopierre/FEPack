numKpoints = 256;
kpars = linspace(-pi, pi, numKpoints);
figure;

% fid1 = fopen('outputs/bulk_spectrum.txt', 'w+');
% fid2 = fopen('outputs/edge_spectrum.txt', 'w+');
% N = 2;
% msg = ['kpar\tminB1\tmaxB1' strjoin(arrayfun(@(x) sprintf('\tminB%d\tmaxB%d', x, x), 2:N, 'UniformOutput', false), '') '\n'];
% fprintf(fid1, msg);
% fprintf(fid2, 'kpar\tEigval');

plot_escurve(kpars);

% for idK = 1:numKpoints
%   kpar = kpars(idK);
% 
%   load(['outputs/out/outputs/EScurves_', num2str(idK)]);
% 
%   [EScurves.XnanBB, EScurves.Xpeaks] = bt.find_nan_boundary_and_peaks(1e5, 1e2);
% 
%   Ycol = EScurves.XnanBB.';
%   Ycol = Ycol(:).';
%   fprintf(fid1, ['%0.5e', repmat('\t%0.5e', 1, numel(Ycol)), '\n'], [kpar, Ycol].');
% 
% 
%   if ~isempty(EScurves.Xpeaks)
% 
%     % plot_escurve(kpars);
%     h = plot(kpar, EScurves.Xpeaks(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
%     fprintf('Save? ''y'' or ''n'': ');
%     choice = input('', 's');
%     if (lower(choice) == 'y')
%       fprintf(fid2, '\n%0.5e\t%0.5e', [kpar, EScurves.Xpeaks(1)]);
%       fprintf('Saved.\n');
%     else
%       fprintf('Skipped.\n');
%     end
%     pause(0.1);
%     delete(h);
% 
%     for idI = 2:numel(EScurves.Xpeaks)
%       % plot_escurve(kpars);
%       h = plot(kpar, EScurves.Xpeaks(idI), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
%       fprintf('Save? ''y'' or ''n'': ');
%       choice = input('', 's');
%       if (lower(choice) == 'y')
%         fprintf(fid2, '\n%0.5e\t%0.5e', [kpar, EScurves.Xpeaks(idI)]);
%         fprintf('Saved.\n');
%       else
%         fprintf('Skipped.\n');
%       end
%       pause(0.1);
%       delete(h);
%     end
%   end
% 
%   for idI = 1:size(EScurves.XnanBB, 1)
%     plot([kpar, kpar], [EScurves.XnanBB(idI, 1), EScurves.XnanBB(idI, 2)], 'Color', [1 1 1]*100/255);
%   end
% 
%   for idI = 1:numel(EScurves.Xpeaks)
%     plot(kpar, EScurves.Xpeaks(idI), 'ro');
%   end
% end
% 
% fclose(fid1);
% fclose(fid2);

function plot_escurve(kpars)
  numKpoints = numel(kpars);

  for idK = 1:numKpoints
    kpar = kpars(idK);

    load(['outputs/out/outputs/EScurves_', num2str(idK)]);
    
    [EScurves.XnanBB, EScurves.Xpeaks] = bt.find_nan_boundary_and_peaks(1e5, 1e2);

    for idI = 1:size(EScurves.XnanBB, 1)
      plot([kpar, kpar], [EScurves.XnanBB(idI, 1), EScurves.XnanBB(idI, 2)], 'Color', [1 1 1]*100/255); 
      hold on;
    end

    for idI = 1:numel(EScurves.Xpeaks)
      plot(kpar, EScurves.Xpeaks(idI), 'ro');
    end
  end
  axis([-pi, pi, -8, 12]);
end