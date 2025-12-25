function construct_solution_3D_par_2(first_use, meshXY, meshLineZ, BCstruct, opts, name_trace)

  if (nargin < 6)
    name_trace = 'sol_trace_Floquet_';
  end

  % Initialization
  sLine = meshLineZ.domain('volumic');
  Ns = meshLineZ.numPoints;
  cutslope = opts.cutvec(2) / opts.cutvec(1);
  FEPackmod = @(x) FEPack.tools.mymod(x, 0, opts.period);

  % Edges
  edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
  edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
  edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
  
  % Basis functions on face X = cst and S/Z-line  
  spBX = BCstruct.spBX;                      % X = cst
  spBS = BCstruct.spBS; NbS = spBS.numBasis; % S/Z-line
  NbY  = NbS * (N0y - 2);

  % Define shear maps
  % /////////////////
  % X = cst
  cut1_0x = opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2);
  cut2_0x = opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2);
  cut1_1x = opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2);
  cut2_1x = opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2);
  %
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([cut1_0x, FEPackmod(cut2_0x + s), zeros(N0x, 1)]);
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([cut1_1x, FEPackmod(cut2_1x + s), zeros(N1x, 1)]);
  %
  % S/Z-line
  [IdSX_S, IdSX_X] = ind2sub([NbS, N0y-2], 1:NbY);
  phisY = eye(N0y-2);
  shearmap.fun0y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPackmod(s           ) 0, 0], IdSX_S);
  shearmap.fun1y = @(s) phisY(:, IdSX_X) .* spBS.evaluateBasisFunctions([FEPackmod(s + cutslope) 0, 0], IdSX_S);
  
  %% Traces of FB transforms
  suffix = opts.suffix;
  numCellsSemiInfinite = opts.numCellsSemiInfinite;
  numCellsInfinite = opts.numCellsInfinite;
  nomdossier = opts.nomdossier;

  if (first_use)
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
  end

  %% Compute inverse FB transform
  idsY = (1:2*numCellsInfinite)';
  sVars = meshLineZ.points(:, 1) + (idsY(:).' - numCellsInfinite - 1) * cutslope;
  sVars = sVars(:);
  
  FloquetPoints = opts.FloquetPoints;
  ptsY = meshXY.points(:, 2);
  cutvec = opts.cutvec;
  period = opts.period;
  
  if (first_use)
    % Locate transverse variable in mesh
    structLoc = sLine.locateInDomain([FEPackmod(sVars), zeros(length(sVars), 2)]);
    elts = sLine.elements(structLoc.elements, :);
    coos = structLoc.barycoos;

    if (opts.numFloquetPoints == 1)
      W = 1;
    else
      W = (FloquetPoints(2) - FloquetPoints(1)) * sqrt(opts.period / (2*pi));
    end

    for idY = 1:2*numCellsInfinite
      
      UsolY = cell(numCellsSemiInfinite, 1);

      for idS = 1:Ns
        fprintf('%d sur %d et %d sur %d\n', idY, 2*numCellsInfinite, idS, Ns);

        I_SY = sub2ind([Ns, 2*numCellsInfinite], idS, idY);
        sVar = sVars(I_SY);
        coosY = coos(I_SY, :);
        eltsY = elts(I_SY, :);
        Usol = zeros(meshXY.numPoints, numCellsSemiInfinite);
        
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

          exp_floquet = exp(1i * FloquetVar * (cutvec(1) * ptsY + (idY - numCellsInfinite - 1) * period)); %#ok
          Usol = Usol + TFB_U .* exp_floquet;
        end

        Usol = Usol * W; %#ok
        
        for idX = 1:numCellsSemiInfinite
          UsolY{idX}(:, idS) = Usol(:, idX);
        end

      end

      for idX = 1:numCellsSemiInfinite
        Usol = UsolY{idX};
        save([nomdossier, 'solution_tmp_', suffix, '_X_', num2str(idX), '_Y_', num2str(idY)], 'Usol', '-v7.3');
      end 

    end
  end

  if ~(first_use)
    % Extract key variables
    X3D = opts.sol3D.points(:, 1);
    Y3D = opts.sol3D.points(:, 2);
    Z3D = opts.sol3D.points(:, 3);
    N3D = length(Z3D);
    idCellsY = opts.sol3D.idCellsY + (numCellsInfinite + 1);        
    idCellsX = abs(opts.sol3D.idCellsX)+1;
    numCellsX = length(idCellsX);
    numCellsY = length(idCellsY);

    % Transverse and directional points and their locations
    ptsTrs = FEPackmod(Z3D - (Y3D + opts.sol3D.idCellsY(:).') * cutslope);
    ptsDir = Y3D / cutvec(1) + zeros(1, numCellsY);
    X3D    = X3D + zeros(1, numCellsY);

    ptsTrs = ptsTrs(:);
    ptsDir = ptsDir(:);
    X3D    = X3D(:);

    structLocTrs = sLine.locateInDomain([ptsTrs, zeros(N3D*numCellsY, 2)]);
    eltsTrs = sLine.elements(structLocTrs.elements, :);
    coosTrs = structLocTrs.barycoos;

    structLocDir = meshXY.domain("volumic").locateInDomain([X3D, ptsDir, zeros(N3D*numCellsY, 1)]);
    eltsDir = meshXY.domain("volumic").elements(structLocDir.elements, :);
    coosDir = structLocDir.barycoos;

    % Reshape elements and barycentric coordinates 
    I3D = cell(3, 2);
    for iTrs = 1:2
      for iDir = 1:3
        I3D{iDir, iTrs} = sub2ind([meshXY.numPoints, Ns], eltsDir(:, iDir), eltsTrs(:, iTrs));
        I3D{iDir, iTrs} = reshape(I3D{iDir, iTrs}, N3D, numCellsY);
      end
    end

    coosTrs = {reshape(coosTrs(:, 1), N3D, numCellsY);...
               reshape(coosTrs(:, 2), N3D, numCellsY)};

    coosDir = {reshape(coosDir(:, 1), N3D, numCellsY);...
               reshape(coosDir(:, 2), N3D, numCellsY);...
               reshape(coosDir(:, 3), N3D, numCellsY)};

    for idY = 1:numCellsY
      for idX = 1:numCellsX
        fprintf('%d sur %d et %d sur %d\n', idY, numCellsY, idX, numCellsX);
        soltmp = load([nomdossier, 'solution_tmp_', suffix,...
                                   '_X_', num2str(idCellsX(idX)),...
                                   '_Y_', num2str(idCellsY(idY))]);
        Usol = zeros(N3D, 1);
        
        for iTrs = 1:2
          for iDir = 1:3
            Usol = Usol + coosDir{iDir}(:, idY) .*...
                          coosTrs{iTrs}(:, idY) .*...
                          soltmp.Usol(I3D{iDir,iTrs}(:,idY));
          end
        end

        save([nomdossier, 'solution_3D_', suffix, '_X_', num2str(idX), '_Y_', num2str(idY)], 'Usol', '-v7.3');
      end
    end 
  end
end

