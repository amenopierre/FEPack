function construct_solution_par(orientation, meshXY, meshLineZ, BCstruct, opts, name_trace)

  if (nargin < 6)
    name_trace = 'sol_trace_Floquet_';
  end

  % Initialization
  sLine = meshLineZ.domain('volumic');
  cutslope = opts.cutvec(2) / opts.cutvec(1);
  spBX = BCstruct.spBX;
  spBY = BCstruct.spBY; 

  % Edges
  edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
  edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
  edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
  edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

  % Define shear map
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2),...
                                                     FEPack.tools.mymod(opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2) + s, 0, opts.period), zeros(N0x, 1)]);
  %
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2),...
                                                     FEPack.tools.mymod(opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2) + s, 0, opts.period), zeros(N1x, 1)]);
  %
  shearmap.fun0y = @(s) spBY.evaluateBasisFunctions([FEPack.tools.mymod(           s * ones(N0y-2, 1), 0, opts.period), meshXY.points(edge0y.IdPoints(2:end-1), 1), zeros(N0y-2, 1)]);
  shearmap.fun1y = @(s) spBY.evaluateBasisFunctions([FEPack.tools.mymod(cutslope + s * ones(N1y-2, 1), 0, opts.period), meshXY.points(edge1y.IdPoints(2:end-1), 1), zeros(N1y-2, 1)]);
  
  if (opts.numFloquetPoints == 1)
    W = 1;
  else
    W = (2*pi/opts.period) ./ (opts.numFloquetPoints - 1);
  end

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
  
  for idY = 1:2*numCellsInfinite
    
    fprintf('%d sur %d\n', idY, 2*numCellsInfinite);
    sVar = sVars(idY);
    Usol = zeros(meshXY.numPoints, opts.numCellsSemiInfinite);
    coosY = coos(idY, :);
    eltsY = elts(idY, :);
   
    for idFB = 1:opts.numFloquetPoints
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
      % r_xi_fun = 1i*FloquetVar*opts.cutvec(1) + 2*pi*(FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2)).';
      % l_xi_fun = sqrt(r_xi_fun.^2 - opts.omega^2);

      % E0 = sinh((1 - orientation*meshXY.points(:, 1)) * l_xi_fun.') .* exp(2i*pi* meshXY.points(:, 2) * (FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2))) .* exp(2i*pi*FourierIdsY(Iy)*sVar) ./ (ones(meshXY.numPoints, 1) * sinh(l_xi_fun.'));
      % %
      % E1 = sinh(orientation*meshXY.points(:, 1) * l_xi_fun.') .* exp(2i*pi * meshXY.points(:, 2) * (FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2))) .* exp(2i*pi*FourierIdsY(Iy)*sVar) ./ (ones(meshXY.numPoints, 1) * sinh(l_xi_fun.'));
      % %
      % TFB_U = E0 * soltrace.R0phi + E1 * soltrace.R1phi;
      % %
      % %
      % % TFB_U = E0 * soltrace.vec + E1 * solguide.R * soltrace.vec;
      % %
      % % TFB_U = (E0 + E1 * exp(-1i*opts.omega)) * soltrace.vec;
      % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      exp_floquet = exp(1i * FloquetVar * (cutvec(1) * ptsY + (idY - numCellsInfinite - 1) * opts.period)); %#ok
      Usol = Usol + TFB_U .* exp_floquet;
    end

    % Usol = Usol * W * sqrt(period / (2*pi)); %#ok

    save([nomdossier, 'solution_', suffix, '_Y_', num2str(idY)], 'Usol', '-v7.3');

  end

  if isfield(opts, 'compute3Dsolution') && (opts.compute3Dsolution)
    % J'utilise le même maillage 2D sur (0, L)² pour construire la trace de la solution
    % sur les faces X = cte, Y = cte, et Z = cte.
    %
    % Sur X = cte, le maillage est orienté suivant YZ
    % Sur Y = cte, le maillage est orienté suivant ZX
    % Sur Z = cte, le maillage est orienté suivant XY
    %
    idsY = (1:2*numCellsInfinite)';
    trsXcst = ones(2*numCellsInfinite, 1) * ...
              FEPack.tools.mymod(opts.meshYZ.points(:, 2)' -...
                                 opts.meshYZ.points(:, 1)'*cutslope, 0, opts.period); % NcellY x Npts
    trsYcst = ones(2*numCellsInfinite, 1) * opts.meshZX.points(:, 1)'; % NcellsY x Npts
    trsZcst = ones(2*numCellsInfinite, 1) * ...
              FEPack.tools.mymod(-opts.meshXY.points(:, 2)'*cutslope, 0, opts.period); % NcellY x Npts

    % Locate transverse variable in mesh
    structLocTrsX = sLine.locateInDomain([trsXcst(:), zeros(numel(trsXcst), 2)]);
    structLocTrsY = sLine.locateInDomain([trsYcst(:), zeros(numel(trsYcst), 2)]);
    structLocTrsZ = sLine.locateInDomain([trsZcst(:), zeros(numel(trsZcst), 2)]);

    eltsTrsX = sLine.elements(structLocTrsX.elements, :); coosTrsX = structLocTrsX.barycoos;
    eltsTrsY = sLine.elements(structLocTrsY.elements, :); coosTrsY = structLocTrsY.barycoos;
    eltsTrsZ = sLine.elements(structLocTrsZ.elements, :); coosTrsZ = structLocTrsZ.barycoos;

    % Locate 3D mesh points in family of 2D meshes
    dirXcst = [zeros(opts.meshYZ.numPoints, 1), opts.meshYZ.points(:, 1)/opts.cutvec(1), zeros(opts.meshYZ.numPoints, 1)]; % Npts x 3
    dirYcst = [opts.meshZX.points(:, 2), zeros(opts.meshZX.numPoints, 2)]; % Npts x 3
    dirZcst = [opts.meshXY.points(:, 1), opts.meshXY.points(:, 2)/opts.cutvec(1), zeros(opts.meshXY.numPoints, 1)]; % Npts x 3

    structLocDirX = meshXY.domain('volumic').locateInDomain(dirXcst);
    structLocDirY = meshXY.domain('volumic').locateInDomain(dirYcst);
    structLocDirZ = meshXY.domain('volumic').locateInDomain(dirZcst);
    
    eltsDirX = meshXY.triangles(structLocDirX.elements, :); coosDirX = structLocDirX.barycoos;
    eltsDirY = meshXY.triangles(structLocDirY.elements, :); coosDirY = structLocDirY.barycoos;
    eltsDirZ = meshXY.triangles(structLocDirZ.elements, :); coosDirZ = structLocDirZ.barycoos;
    
    FloquetPoints = opts.FloquetPoints;
    numFloquetPoints = opts.numFloquetPoints;
    f_X = opts.meshYZ.points(:, 1)' + (idsY - numCellsInfinite - 1);
    f_Y = (idsY - numCellsInfinite - 1) * ones(1, opts.meshZX.numPoints);
    f_Z = opts.meshXY.points(:, 2)' + (idsY - numCellsInfinite - 1);
    
    period = opts.period;
    numCellsSemiInfinite = opts.numCellsSemiInfinite;
    % idCpt = 0;
    
    for idI = 1:2*numCellsInfinite 
      fprintf('%d sur %d\n', idI, 2*numCellsInfinite);

      %% X = cte
      %  ///////
      U3DsolX = zeros(opts.meshYZ.numPoints, numCellsSemiInfinite);
      NX = opts.meshYZ.numPoints;

      id_IJ = sub2ind([2*numCellsInfinite, NX], idI*ones(1, NX), 1:NX);
      coosTrs_I = coosTrsX(id_IJ, :);
      eltsTrs_I = eltsTrsX(id_IJ, :);
      f_X_I = f_X(idI, :);

      parfor idJ = 1:NX
        % idCpt = idCpt + 1;
        % fprintf('%d sur %d\n', idCpt, 2*numCellsInfinite*(opts.meshYZ.numPoints+opts.meshZX.numPoints+opts.meshXY.numPoints));
        
        coosDir_J = coosDirX(idJ, :); 
        eltsDir_J = eltsDirX(idJ, :);
        trsVar = trsXcst(idI, idJ);

        coosTrs_IJ = coosTrs_I(idJ, :); 
        eltsTrs_IJ = eltsTrs_I(idJ, :);
    
        f_IJ = f_X_I(idJ);
        
        Usol = zeros(1, numCellsSemiInfinite);
        
        for idFB = 1:numFloquetPoints
          FloquetVar = FloquetPoints(idFB); %#ok
          solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
          soltrace = load([nomdossier, name_trace, num2str(idFB)]);
          
          e0x = 0; e1x = 0; e0y = 0; e1y = 0; 
          for iTrs = 1:2
            sollcell = load([nomdossier, 'local_cell_sol_', suffix,...
                            '_Floquet_', num2str(idFB), '_S_', num2str(eltsTrs_IJ(iTrs))]); % #ok
            for iDir = 1:3
              e0x = e0x + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E0x(eltsDir_J(iDir), :); % #ok
              e1x = e1x + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E1x(eltsDir_J(iDir), :);
              e0y = e0y + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E0y(eltsDir_J(iDir), :);
              e1y = e1y + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E1y(eltsDir_J(iDir), :);
            end
          end

          TFB_U = e0x * shearmap.fun0x(trsVar) * soltrace.R0phi... % X direction, E0
                + e1x * shearmap.fun1x(trsVar) * soltrace.R1phi... % X direction, E1
                + e0y * shearmap.fun0y(trsVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
                + e1y * shearmap.fun1y(trsVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
                + e0y * shearmap.fun0y(trsVar) * solguide.DtD1 * soltrace.R1phi...    % Y direction, E1
                + e1y * shearmap.fun1y(trsVar) * solguide.DtD1 * soltrace.R1phi; %#ok % Y direction, E1

          Usol = Usol + TFB_U .* exp(1i * FloquetVar * f_IJ * period);
        end

        U3DsolX(idJ, :) = Usol * W * sqrt(period / (2*pi));
      end
      Usol = U3DsolX; %#ok 
      save([nomdossier, 'solution3D_Xcst_', suffix, '_Y_', num2str(idI)], 'Usol', '-v7.3');

      %% Y = cte
      %  ///////
      U3DsolY = zeros(opts.meshZX.numPoints, numCellsSemiInfinite);
      NY = opts.meshZX.numPoints;

      id_IJ = sub2ind([2*numCellsInfinite, NY], idI*ones(1, NY), 1:NY);
      coosTrs_I = coosTrsY(id_IJ, :);
      eltsTrs_I = eltsTrsY(id_IJ, :);
      f_Y_I = f_Y(idI, :);

      parfor idJ = 1:NY
        % idCpt = idCpt + 1;
        % fprintf('%d sur %d\n', idCpt, 2*numCellsInfinite*(opts.meshZX.numPoints+opts.meshZX.numPoints+opts.meshXY.numPoints));
        
        coosDir_J = coosDirY(idJ, :); 
        eltsDir_J = eltsDirY(idJ, :);
        trsVar = trsYcst(idI, idJ);

        coosTrs_IJ = coosTrs_I(idJ, :); 
        eltsTrs_IJ = eltsTrs_I(idJ, :);
    
        f_IJ = f_Y_I(idJ);
        
        Usol = zeros(1, numCellsSemiInfinite);
        
        for idFB = 1:numFloquetPoints
          FloquetVar = FloquetPoints(idFB); %#ok
          solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
          soltrace = load([nomdossier, name_trace, num2str(idFB)]);
          
          e0x = 0; e1x = 0; e0y = 0; e1y = 0; 
          for iTrs = 1:2
            sollcell = load([nomdossier, 'local_cell_sol_', suffix,...
                            '_Floquet_', num2str(idFB), '_S_', num2str(eltsTrs_IJ(iTrs))]); % #ok
            for iDir = 1:3
              e0x = e0x + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E0x(eltsDir_J(iDir), :); % #ok
              e1x = e1x + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E1x(eltsDir_J(iDir), :);
              e0y = e0y + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E0y(eltsDir_J(iDir), :);
              e1y = e1y + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E1y(eltsDir_J(iDir), :);
            end
          end

          TFB_U = e0x * shearmap.fun0x(trsVar) * soltrace.R0phi... % X direction, E0
                + e1x * shearmap.fun1x(trsVar) * soltrace.R1phi... % X direction, E1
                + e0y * shearmap.fun0y(trsVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
                + e1y * shearmap.fun1y(trsVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
                + e0y * shearmap.fun0y(trsVar) * solguide.DtD1 * soltrace.R1phi...    % Y direction, E1
                + e1y * shearmap.fun1y(trsVar) * solguide.DtD1 * soltrace.R1phi; %#ok % Y direction, E1

          Usol = Usol + TFB_U .* exp(1i * FloquetVar * f_IJ * period);
        end

        U3DsolY(idJ, :) = Usol * W * sqrt(period / (2*pi));
      end
      Usol = U3DsolY; %#ok 
      save([nomdossier, 'solution3D_Ycst_', suffix, '_Y_', num2str(idI)], 'Usol', '-v7.3');

      %% Z = cte
      %  ///////
      U3DsolZ = zeros(opts.meshXY.numPoints, numCellsSemiInfinite);
      NZ = opts.meshXY.numPoints;

      id_IJ = sub2ind([2*numCellsInfinite, NZ], idI*ones(1, NZ), 1:NZ);
      coosTrs_I = coosTrsZ(id_IJ, :);
      eltsTrs_I = eltsTrsZ(id_IJ, :);
      f_Z_I = f_Z(idI, :);

      parfor idJ = 1:NZ
        % idCpt = idCpt + 1;
        % fprintf('%d sur %d\n', idCpt, 2*numCellsInfinite*(opts.meshXY.numPoints+opts.meshXY.numPoints+opts.meshXY.numPoints));
        
        coosDir_J = coosDirZ(idJ, :); 
        eltsDir_J = eltsDirZ(idJ, :);
        trsVar = trsZcst(idI, idJ);

        coosTrs_IJ = coosTrs_I(idJ, :); 
        eltsTrs_IJ = eltsTrs_I(idJ, :);
    
        f_IJ = f_Z_I(idJ);
        
        Usol = zeros(1, numCellsSemiInfinite);
        
        for idFB = 1:numFloquetPoints
          FloquetVar = FloquetPoints(idFB); %#ok
          solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
          soltrace = load([nomdossier, name_trace, num2str(idFB)]);
          
          e0x = 0; e1x = 0; e0y = 0; e1y = 0; 
          for iTrs = 1:2
            sollcell = load([nomdossier, 'local_cell_sol_', suffix,...
                            '_Floquet_', num2str(idFB), '_S_', num2str(eltsTrs_IJ(iTrs))]); % #ok
            for iDir = 1:3
              e0x = e0x + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E0x(eltsDir_J(iDir), :); % #ok
              e1x = e1x + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E1x(eltsDir_J(iDir), :);
              e0y = e0y + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E0y(eltsDir_J(iDir), :);
              e1y = e1y + coosTrs_IJ(iTrs) * coosDir_J(iDir) * sollcell.E1y(eltsDir_J(iDir), :);
            end
          end

          TFB_U = e0x * shearmap.fun0x(trsVar) * soltrace.R0phi... % X direction, E0
                + e1x * shearmap.fun1x(trsVar) * soltrace.R1phi... % X direction, E1
                + e0y * shearmap.fun0y(trsVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
                + e1y * shearmap.fun1y(trsVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
                + e0y * shearmap.fun0y(trsVar) * solguide.DtD1 * soltrace.R1phi...    % Y direction, E1
                + e1y * shearmap.fun1y(trsVar) * solguide.DtD1 * soltrace.R1phi; %#ok % Y direction, E1

          Usol = Usol + TFB_U .* exp(1i * FloquetVar * f_IJ * period);
        end

        U3DsolZ(idJ, :) = Usol * W * sqrt(period / (2*pi));
      end
      Usol = U3DsolZ; %#ok 
      save([nomdossier, 'solution3D_Zcst_', suffix, '_Y_', num2str(idI)], 'Usol', '-v7.3');
    end
  end
end

