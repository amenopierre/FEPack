function construct_solution_par(orientation, meshXY, meshLineZ, BCstruct, opts, name_trace)

  if (nargin < 6)
    name_trace = 'sol_trace_Floquet_';
  end

  % Initialization
  sLine = meshLineZ.domain('volumic');
  cutslope = opts.cutvec(2) / opts.cutvec(1);
  
  % Edges
  edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
  edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
  edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
  edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

  % Basis functions on face X = cst and S/Z-line  
  spBX = BCstruct.spBX; NbX = spBX.numBasis; % X = cst
  spBS = BCstruct.spBS; NbS = spBS.numBasis; % S/Z-line
  NbY = NbS * (N0y - 2);

  % Define shear maps
  % /////////////////
  % X = cst
  cut1_0x = opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2);
  cut2_0x = opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2);
  cut1_1x = opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2);
  cut2_1x = opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2);
  %
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([cut1_0x, FEPack.tools.mymod(cut2_0x + s, 0, opts.period), zeros(N0x, 1)]);
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([cut1_1x, FEPack.tools.mymod(cut2_1x + s, 0, opts.period), zeros(N1x, 1)]);
  %
  % S/Z-line
  [IdSX_S, IdSX_X] = ind2sub([NbS, N0y-2], 1:NbY);
  phisY = eye(N0y-2);
  shearmap.fun0y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPack.tools.mymod(s,            0, opts.period) 0, 0], IdSX_S);
  shearmap.fun1y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPack.tools.mymod(s + cutslope, 0, opts.period) 0, 0], IdSX_S);
  
  %% Traces of FB transforms
  suffix = opts.suffix;
  numCellsSemiInfinite = opts.numCellsSemiInfinite;
  numCellsInfinite = opts.numCellsInfinite;
  nomdossier = opts.nomdossier;

  parfor idFB = 1:opts.numFloquetPoints
    soltrace = load([nomdossier, name_trace, num2str(idFB)]);
    solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
    
    soltrace.R0phi = zeros(size(soltrace.vec, 1), numCellsSemiInfinite);
    soltrace.R0phi(:, 1) = soltrace.vec;

    for idX = 2:numCellsSemiInfinite
      soltrace.R0phi(:, idX) = solguide.R * soltrace.R0phi(:, idX-1);
    end
    soltrace.R1phi = solguide.D * soltrace.R0phi;

    % save([nomdossier, name_trace, num2str(idFB)], 'R0phi', 'R1phi', '-append');
    FEPack.tools.parsave([nomdossier, name_trace, num2str(idFB)], soltrace, true);
  end

  %% Compute inverse FB transform
  idsY = (1:2*numCellsInfinite)';
  sVars = (idsY - numCellsInfinite - 1) * cutslope;  % Transverse variables
  
  % Locate transverse variable in mesh
  structLoc = sLine.locateInDomain([FEPack.tools.mymod(sVars, 0, opts.period), zeros(2*numCellsInfinite, 2)]);
  elts = sLine.elements(structLoc.elements, :);
  coos = structLoc.barycoos;

  FloquetPoints = opts.FloquetPoints;
  ptsY = meshXY.points(:, 2);
  cutvec = opts.cutvec;
  period = opts.period;
  
  if (opts.numFloquetPoints == 1)
    W = 1;
  else
    W = (FloquetPoints(2) - FloquetPoints(1)) * sqrt(opts.period / (2*pi));
  end

  fid = fopen([nomdossier, 'trace_solution_X_0_', suffix, '.txt'], 'w+');
  trace_sol = [];

  for idY = 1:2*numCellsInfinite
    
    fprintf('%d sur %d\n', idY, 2*numCellsInfinite);
    sVar = sVars(idY);
    Usol = zeros(meshXY.numPoints, opts.numCellsSemiInfinite);
    coosY = coos(idY, :);
    eltsY = elts(idY, :);
   
    parfor idFB = 1:opts.numFloquetPoints
      FloquetVar = FloquetPoints(idFB);
      solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
      soltrace = load([nomdossier, name_trace, num2str(idFB)]);
      sc1 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(eltsY(1))]); %#ok
      sc2 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(eltsY(2))]);
      
      e0x = coosY(1) * sc1.E0x + coosY(2) * sc2.E0x; %#ok
      e1x = coosY(1) * sc1.E1x + coosY(2) * sc2.E1x;
      e0y = coosY(1) * sc1.E0y + coosY(2) * sc2.E0y;
      e1y = coosY(1) * sc1.E1y + coosY(2) * sc2.E1y;
      
      TFB_U = e0x * shearmap.fun0x(sVar) * soltrace.R0phi... % X direction, E0
            + e1x * shearmap.fun1x(sVar) * soltrace.R1phi... % X direction, E1
            + e0y * shearmap.fun0y(sVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
            + e1y * shearmap.fun1y(sVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
            + e0y * shearmap.fun0y(sVar) * solguide.DtD1 * soltrace.R1phi...    % Y direction, E1
            + e1y * shearmap.fun1y(sVar) * solguide.DtD1 * soltrace.R1phi; %#ok % Y direction, E1
      
      % % Cas homogène
      % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % FourierIdsX = spBX.FourierIds.X(:).'; dx = length(FourierIdsX);
      % FourierIdsY = spBX.FourierIds.Y(:).'; dy = length(FourierIdsY);
      % FourierIdsZ = spBX.FourierIds.Z(:).'; dz = length(FourierIdsZ);

      % [Ix, Iy, ~] = ind2sub([dx, dy, dz], (1:spBX.numBasis)');
      % r_xi_fun = FloquetVar*opts.cutvec(1) + 2*pi*(FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2)).';
      % l_xi_fun = sqrt(r_xi_fun.^2 - opts.omega^2);

      % E0 = sinh((1 - orientation*meshXY.points(:, 1)) * l_xi_fun.') .* exp(2i*pi* meshXY.points(:, 2) * (FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2))) .* exp(2i*pi*FourierIdsY(Iy)*sVar) ./ (ones(meshXY.numPoints, 1) * sinh(l_xi_fun.'));
      % %
      % E1 = sinh(orientation*meshXY.points(:, 1) * l_xi_fun.') .* exp(2i*pi * meshXY.points(:, 2) * (FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2))) .* exp(2i*pi*FourierIdsY(Iy)*sVar) ./ (ones(meshXY.numPoints, 1) * sinh(l_xi_fun.'));
      % %
      % %
      % TFB_Uvec = exp(-orientation * meshXY.points(:, 1) * l_xi_fun.') .* ...
      %            exp(2i*pi*meshXY.points(:, 2)*(FourierIdsX(Ix)*opts.cutvec(1) + FourierIdsY(Iy)*opts.cutvec(2))) .*...
      %            exp(2i*pi*FourierIdsY(Iy)*sVar);

      % TFB_U = E0 * soltrace.R0phi + E1 * soltrace.R1phi;
      % % TFB_U = E0 * soltrace.vec + E1 * solguide.R * soltrace.vec;
      % % TFB_U = E0 * soltrace.vec + E1 * diag(exp(-l_xi_fun)) * soltrace.vec;
      % % TFB_U = TFB_Uvec * soltrace.vec;
      % %
      % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      exp_floquet = exp(1i * FloquetVar * (cutvec(1) * ptsY + (idY - numCellsInfinite - 1) * period)); %#ok
      Usol = Usol + TFB_U .* exp_floquet;
    end

    Usol = Usol * W; %#ok

    save([nomdossier, 'solution_', suffix, '_Y_', num2str(idY)], 'Usol', '-v7.3');

    % Save trace of solution at interface (x = 0)
    trace_sol = [trace_sol; [ptsY(edge0x.IdPoints) + (idY - numCellsInfinite - 1) * period / opts.cutvec(1), real(Usol(edge0x.IdPoints, 1)), imag(Usol(edge0x.IdPoints, 1))]];

  end

  [~, Ikey] = unique(trace_sol(:, 1));
  trace_sol = trace_sol(Ikey, :);
  fprintf(fid, '%0.5e\t%0.5e\t%0.5e\n', trace_sol.');
  fclose(fid);

end