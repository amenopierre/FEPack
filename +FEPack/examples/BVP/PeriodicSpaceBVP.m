function U = PeriodicSpaceBVP(semiInfiniteDirection, infiniteDirections,...
                              volBilinearIntg_pos, mesh_pos, BCstruct_pos, numCellsSemiInfinite_pos,...
                              volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCellsSemiInfinite_neg,...
                              volBilinearIntg_int, volLinearIntg, mesh_int, spBint_pos, spBint_neg,...
                              numCellsInfinite, numFloquetPoints, opts)
  % PeriodicSpaceBVP

  %% % ************************* %
  %  % Preliminary verifications %
  %  % ************************* %
  pgbars = findall(0,'type','figure','tag','TMWWaitbar');
  delete(pgbars);

  infiniteDirections = unique(infiniteDirections(:).');
  opts.verbose = 0;

  if min((1:mesh_pos.dimension) ~= semiInfiniteDirection)
    % The component along which the domain is semi-infinite should be between
    % 1 and the mesh dimension
    error('Les sous-domaines ne peuvent pas être semi-infinis dans la direction %d.', semiInfiniteDirection);
  end

  for idI = 1:length(infiniteDirections)
    if min((1:mesh_pos.dimension) ~= infiniteDirections(idI))
      % The components along which the domain is infinite should be between
      % 1 and the mesh dimension
      error('Le domaine ne peut pas être infini dans la direction %d.', infiniteDirection(idI));
    end
  end

  if max(infiniteDirections == semiInfiniteDirection)
    % The half-space cannot be both infinite and semi-infinite in the same direction
    error('Le domaine ne peut pas être à la fois infini et semi-infini dans la direction %d.', semiInfiniteDirection)
  end

  if (length(semiInfiniteDirection) > 1)
    % The space cannot be semi-infinte in more than one direction
    error('Le domaine ne peut pas être semi-infini dans plus d''une direction.');
  end

  if (length(numCellsInfinite) ~= length(infiniteDirections))
    % The length of the numCells should be equal to the dimension
    error('Le vecteur numCellsInfinite doit avoir %d composantes.', length(infiniteDirections));
  end

  if (length(numFloquetPoints) ~= length(infiniteDirections))
    % The length of the numFloquetPoints vector should be equal to the dimension
    error('Le vecteur numFloquetPoints doit avoir %d composantes.', length(infiniteDirections));
  end

  if max(numFloquetPoints < 1)
    error('Il faut au moins 2 points pour la variable de Floquet.');
  end

  %% % ********************************** %
  %  % Compute the elementary FE matrices %
  %  % ********************************** %
  Ni = length(infiniteDirections);

  % Positive direction
  Npos = mesh_pos.numPoints;
  AAelem_pos = cell(Ni +1, Ni + 1);
  for idU = 1:(Ni + 1)
    for idV = 1:(Ni + 1)
      AAelem_pos{idV, idU} = sparse(Npos, Npos);
    end
  end

  for idT = 1:length(volBilinearIntg_pos.alpha_u)
    fun = volBilinearIntg_pos.fun{idT};

    Au = cell(Ni + 1, 1);
    Av = cell(Ni + 1, 1);
    Au{1} = volBilinearIntg_pos.alpha_u{idT};
    Av{1} = volBilinearIntg_pos.alpha_v{idT};
    for idB = 1:Ni
      Au{idB + 1} = Au{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
      Av{idB + 1} = Av{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
    end

    for idU = 1:(Ni + 1)
      for idV = 1:(Ni + 1)
        AAelem_pos{idV, idU} = AAelem_pos{idV, idU} + FEPack.pdes.Form.global_matrix(mesh_pos.domain('volumic'), Au{idU}, Av{idV}, fun);
      end
    end
  end

  % Negative direction
  Nneg = mesh_neg.numPoints;
  AAelem_neg = cell(Ni +1, Ni + 1);
  for idU = 1:(Ni + 1)
    for idV = 1:(Ni + 1)
      AAelem_neg{idV, idU} = sparse(Nneg, Nneg);
    end
  end

  for idT = 1:length(volBilinearIntg_neg.alpha_u)
    fun = volBilinearIntg_neg.fun{idT};

    Au = cell(Ni + 1, 1);
    Av = cell(Ni + 1, 1);
    Au{1} = volBilinearIntg_neg.alpha_u{idT};
    Av{1} = volBilinearIntg_neg.alpha_v{idT};
    for idB = 1:Ni
      Au{idB + 1} = Au{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
      Av{idB + 1} = Av{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
    end

    for idU = 1:(Ni + 1)
      for idV = 1:(Ni + 1)
        AAelem_neg{idV, idU} = AAelem_neg{idV, idU} + FEPack.pdes.Form.global_matrix(mesh_neg.domain('volumic'), Au{idU}, Av{idV}, fun);
      end
    end
  end

  % Interior domain
  Nint = mesh_int.numPoints;
  AAelem_int = cell(Ni +1, Ni + 1);
  for idU = 1:(Ni + 1)
    for idV = 1:(Ni + 1)
      AAelem_int{idV, idU} = sparse(Nint, Nint);
    end
  end

  for idT = 1:length(volBilinearIntg_int.alpha_u)
    fun = volBilinearIntg_int.fun{idT};

    Au = cell(Ni + 1, 1);
    Av = cell(Ni + 1, 1);
    Au{1} = volBilinearIntg_int.alpha_u{idT};
    Av{1} = volBilinearIntg_int.alpha_v{idT};
    for idB = 1:Ni
      Au{idB + 1} = Au{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
      Av{idB + 1} = Av{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
    end

    for idU = 1:(Ni + 1)
      for idV = 1:(Ni + 1)
        AAelem_int{idV, idU} = AAelem_int{idV, idU} + FEPack.pdes.Form.global_matrix(mesh_int.domain('volumic'), Au{idU}, Av{idV}, fun);
      end
    end
  end

  %% % *************************************************** %
  %  % Compute the Floquet-Bloch transform of the solution %
  %  % *************************************************** %
  numFBpoints = prod(numFloquetPoints);
  FBpoints = cell(Ni, 1);

  % Discretization grid for the Floquet variable
  for idI = 1:Ni
    FBpoints{idI} = linspace(-pi, pi, numFloquetPoints(idI));
  end

  % Linearized to subindices for Floquet variable
  FBdim = [1 1 1];
  FBdim(infiniteDirections) = numFloquetPoints;
  [I1, I2, I3] = ind2sub(FBdim, 1:numFBpoints);
  FloquetIds = [I1; I2; I3];
  FloquetIds = FloquetIds(infiniteDirections, :); % Ni-by-numFBpoints
  TFBU = cell(numFBpoints, 1);

  K = zeros(numFBpoints, Ni);
  for idI = 1:Ni
    K(:, idI) = FBpoints{idI}(FloquetIds(idI, :)).';
  end

  volLinearIntg_FB = copy(volLinearIntg);

  % Solve a family of half-guide problems
  % /////////////////////////////////////
  pgbar = waitbar(0,'','Name','Problèmes de demi-guide',...
                       'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
  for idFB = 1:numFBpoints
    % Check for clicked Cancel button
    if getappdata(pgbar,'canceling')
      delete(pgbar);
      break
    end

    % Update waitbar and message
    waitbar(idFB/numFBpoints, pgbar, 'Problèmes de demi-guide');

    % The FE matrix is a linear combination of the elementary pieces
    FBvar = [1, 1i*K(idFB, :)];
    AApos = sparse(Npos, Npos);
    AAneg = sparse(Nneg, Nneg);
    AAint = sparse(Nint, Nint);
    for idU = 1:(Ni + 1)
      for idV = 1:(Ni + 1)
        AApos = AApos + FBvar(idV)' * AAelem_pos{idV, idU} * FBvar(idU);
        AAneg = AAneg + FBvar(idV)' * AAelem_neg{idV, idU} * FBvar(idU);
        AAint = AAint + FBvar(idV)' * AAelem_int{idV, idU} * FBvar(idU);
      end
    end

    % The Floquet-Bloch transform of the boundary data along the relevant directions
    volLinearIntg_FB.fun = @(x) FEPack.tools.BlochTransform(x, K(idFB, :), volLinearIntg.fun, infiniteDirections);

    % Compute the solution
    TFBU{idFB} = PeriodicGuideBVP(semiInfiniteDirection,...
                                  AApos, mesh_pos, BCstruct_pos, numCellsSemiInfinite_pos,...
                                  AAneg, mesh_neg, BCstruct_neg, numCellsSemiInfinite_neg,...
                                  AAint, volLinearIntg, mesh_int,...
                                  spBint_pos, spBint_neg, opts);
  end

  delete(pgbar);

  %% % ***************************************** %
  %  % Apply the inverse Floquet-Bloch transform %
  %  % ***************************************** %
  % Positive side
  % /////////////
  numCells = [1 1 1];
  numCells(infiniteDirections) = 2*numCellsInfinite;
  numCells(semiInfiniteDirection) = numCellsSemiInfinite_pos;
  Nu = prod(numCells);

  [I1, I2, I3] = ind2sub(numCells, 1:Nu);
  pointsIds = [I1; I2; I3]; % 3-by-Nu
  tau = pointsIds(infiniteDirections, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
  W = prod(2*pi ./ (numFloquetPoints - 1)); % 1-by-1
  U.positive = zeros(Npos, Nu);

  for idFB = 1:numFBpoints
    % The integral that defines the inverse Floquet-Bloch transform is computed
    % using a rectangular rule.
    exp_k_dot_x = exp(1i * mesh_pos.points(:, infiniteDirections) * K(idFB, :).'); % N-by-1
    exp_k_dot_tau = exp(1i * K(idFB, :) * tau); % 1-by-Nu
    U_TFB = TFBU{idFB}.positive(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

    U.positive = U.positive + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
  end

  U.positive = U.positive / sqrt(2*pi);

  % Negative side
  % /////////////
  numCells = [1 1 1];
  numCells(infiniteDirections) = 2*numCellsInfinite;
  numCells(semiInfiniteDirection) = numCellsSemiInfinite_neg;
  Nu = prod(numCells);

  [I1, I2, I3] = ind2sub(numCells, 1:Nu);
  pointsIds = [I1; I2; I3]; % 3-by-Nu
  tau = pointsIds(infiniteDirections, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
  W = prod(2*pi ./ (numFloquetPoints - 1)); % 1-by-1
  U.negative = zeros(Nneg, Nu);

  for idFB = 1:numFBpoints
    % The integral that defines the inverse Floquet-Bloch transform is computed
    % using a rectangular rule.
    exp_k_dot_x = exp(1i * mesh_neg.points(:, infiniteDirections) * K(idFB, :).'); % N-by-1
    exp_k_dot_tau = exp(1i * K(idFB, :) * tau); % 1-by-Nu
    U_TFB = TFBU{idFB}.negative(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

    U.negative = U.negative + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
  end

  U.negative = U.negative / sqrt(2*pi);

  % Interior domain
  % ///////////////
  numCells = [1 1 1]; % Number of cells along each component
  numCells(infiniteDirections) = 2*numCellsInfinite;
  Nu = prod(numCells);

  [I1, I2, I3] = ind2sub(numCells, 1:Nu);
  pointsIds = [I1; I2; I3]; % 3-by-Nu
  tau = pointsIds(infiniteDirections, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
  W = prod(2*pi ./ (numFloquetPoints - 1)); % 1-by-1
  U.interior = zeros(Nint, Nu);

  for idFB = 1:numFBpoints
    % The integral that defines the inverse Floquet-Bloch transform is computed
    % using a rectangular rule.
    exp_k_dot_x = exp(1i * mesh_int.points(:, infiniteDirections) * K(idFB, :).'); % N-by-1
    exp_k_dot_tau = exp(1i * K(idFB, :) * tau); % 1-by-Nu
    U_TFB = TFBU{idFB}.interior(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

    U.interior = U.interior + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
  end

  U.interior = U.interior / sqrt(2*pi);
end
