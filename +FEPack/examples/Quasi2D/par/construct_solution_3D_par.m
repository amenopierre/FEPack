function construct_solution_3D_par(orientation, meshXY, meshLineZ, BCstruct, opts, name_trace)

  % INPUTS: * sol3D.points: N x 3 vector [x y z] with |x|, y, z in (0, 1)
  %         * sol3D.idCellsX: vector of nonnegative integers, that indicates
  %                           the cell indices (|x| ± (n, n+1)) where the 
  %                           solution is computed
  %         * sol3D.idCellsY: vector of integers, indicating the cell indices
  %                           (y + (l, l+1)) where the solution is computed
  
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
  NbY  = NbS * (N0y - 2);

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
  sVars = FEPack.tools.mymod((idsY - numCellsInfinite - 1) * cutslope, 0, opts.period);  % Transverse variables
  
  % Locate transverse variable in mesh
  structLoc = sLine.locateInDomain([sVars, zeros(2*numCellsInfinite, 2)]);
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

  % Extract key variables
  X3D = opts.sol3D.points(:, 1);
  Y3D = opts.sol3D.points(:, 2);
  Z3D = opts.sol3D.points(:, 3);
  N3D = length(Z3D);
  idCellsY = opts.sol3D.idCellsY; numCellsY = length(idCellsY);
  idCellsX = abs(opts.sol3D.idCellsX)+1; numCellsX = length(idCellsX);

  % Transverse and directional points and their locations
  ptsTrs = FEPack.tools.mymod(Z3D - Y3D * cutslope, 0, period);
  ptsDir = Y3D / cutvec(1);
  
  structLocTrs = sLine.locateInDomain([ptsTrs, zeros(N3D, 2)]);
  eltsTrs = sLine.elements(structLocTrs.elements, :);
  coosTrs = structLocTrs.barycoos;

  structLocDir = meshXY.domain("volumic").locateInDomain([X3D, ptsDir, zeros(N3D, 1)]);
  eltsDir = meshXY.domain("volumic").elements(structLocDir.elements, :);
  coosDir = structLocDir.barycoos;

  for idY = 1:numCellsY
    disp(idY)
    U3Dsol = zeros(N3D, numCellsX);

    parfor idI = 1:N3D
      disp(idI)
      eltsTrs_I = eltsTrs(idI, :);
      coosTrs_I = coosTrs(idI, :);
      eltsDir_I = eltsDir(idI, :);
      coosDir_I = coosDir(idI, :);

      Y3Dvar = Y3D(idI);
      trsVar = ptsTrs(idI);
      Usol = zeros(1, numCellsX);
      
      for idFB = 1:opts.numFloquetPoints
        FloquetVar = FloquetPoints(idFB); %#ok
        solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
        soltrace = load([nomdossier, name_trace, num2str(idFB)]);
        
        e0x = 0; e1x = 0; e0y = 0; e1y = 0;

        for iTrs = 1:2
          sollcell = load([nomdossier, 'local_cell_sol_', suffix,...
                          '_Floquet_', num2str(idFB), '_S_', num2str(eltsTrs_I(iTrs))]); % #ok
          for iDir = 1:3
            e0x = e0x + coosTrs_I(iTrs) * coosDir_I(iDir) * sollcell.E0x(eltsDir_I(iDir), :); % #ok
            e1x = e1x + coosTrs_I(iTrs) * coosDir_I(iDir) * sollcell.E1x(eltsDir_I(iDir), :);
            e0y = e0y + coosTrs_I(iTrs) * coosDir_I(iDir) * sollcell.E0y(eltsDir_I(iDir), :);
            e1y = e1y + coosTrs_I(iTrs) * coosDir_I(iDir) * sollcell.E1y(eltsDir_I(iDir), :);
          end
        end

        TFB_U = e0x * shearmap.fun0x(trsVar) * soltrace.R0phi(:, idCellsX)... % X direction, E0
              + e1x * shearmap.fun1x(trsVar) * soltrace.R1phi(:, idCellsX)... % X direction, E1
              + e0y * shearmap.fun0y(trsVar) * solguide.DtD0 * soltrace.R0phi(:, idCellsX)...    % Y direction, E0
              + e1y * shearmap.fun1y(trsVar) * solguide.DtD0 * soltrace.R0phi(:, idCellsX)...    % Y direction, E0
              + e0y * shearmap.fun0y(trsVar) * solguide.DtD1 * soltrace.R1phi(:, idCellsX)...    % Y direction, E1
              + e1y * shearmap.fun1y(trsVar) * solguide.DtD1 * soltrace.R1phi(:, idCellsX); %#ok % Y direction, E1

        Usol = Usol + TFB_U * exp(1i * FloquetVar * (Y3Dvar + idCellsY(idY)) * period);
      end

      U3Dsol(idI, :) = Usol * W;
    end

    for idX = 1:numCellsX
      Usol = U3Dsol(:, idX);
      save([nomdossier, 'solution_3D_', suffix, '_X_', num2str(idX), '_Y_', num2str(idY)], 'Usol', '-v7.3');
    end
  end

end

