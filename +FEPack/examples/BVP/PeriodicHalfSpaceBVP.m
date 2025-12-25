function U = PeriodicHalfSpaceBVP(mesh, orientation, semiInfiniteDirection,...
                                  infiniteDirections, volBilinearIntg,...
                                  BCstruct, numCellsSemiInfinite, numCellsInfinite,...
                                  numFloquetPoints, opts)
  % PeriodicHalfSpaceBVP(mesh, orientation, infiniteDirection, volBilinearIntg, BCstruct, numCells, opts) solves
  %
  %       Find u in V such that tau(u) = g, and
  %                A(u, v) + S(u, v) = Sf(v) for any v in V such that tau(v) = 0,
  %
  % where tau is a 0-order linear operator, in a half-space.
  %
  % INPUTS: * mesh, Mesh object.
  %
  %         * orientation, a number in {-1, 1} that represent whether the domain
  %           is infinite along a positive or a negative axis.
  %
  %         * semiInfiniteDirection, a number that specifies the direction
  %           along which the domain is semi-infinite.
  %
  %         * infiniteDirections, a Ni-sized vector that contains the
  %           directions along which the domain is infinite.
  %
  %         * bilinearIntg, integrand associated to the volumic bilinear form.
  %            For instance, bilinearIntg = grad(u)*grad(v) - (omega^2)*id(u)*id(v).
  %
  %         * BCstruct, a struct with multiple fields:
  %             - spB0, SpectralBasis object, set of basis functions
  %               associated to {x_i = 0} (if the guide is infinite along x_i).
  %             - spB1, same as spB0, but with {x_i = 1}. NOTE that
  %               spB0 and spB1 should correspond to the same
  %               basis on different domains.
  %             - BCdu, BCu and phi are such that the solution satisfies the
  %               boundary condition
  %                            BCdu * (du/dn) + BCu * u = phi.
  %               -> BCdu is a scalar.
  %               -> BCu can be either
  %                   ** a scalar,
  %                   ** a function_handle, or
  %                   ** an Nb-by-Nb matrix which represents an operator in the
  %                      spectral space spanned by the Nb spB0 functions.
  %                      NOTE: in this case, one has to specify if BCu is a
  %                      weak L2 representation or a matrix representation of
  %                      the operator, by means of BCstruct.representation.
  %               -> phi is a function_handle.
  %             - representation: required if BCu is a matrix. It can be set
  %               either to 'weak evaluation' or 'projection' depending on
  %               how BCu is defined:
  %               -> if BCu_{i, j} = <Op(phi_j), phi_i>_L2, then 'weak-evaluation'
  %               -> If Op(phi_j) = BCu_{1, j} phi_1 + ... + BCu_{N, j} phi_N,
  %                  then 'projection'.
  %
  %         * numCellsSemiInfinite, the number of cells in the semi-infinite direction.
  %
  %         * numCellsInfinite, Ni-sized vector containing the number of cells
  %           in the infinite directions.
  %
  %         * numFloquetPoints, Ni-sized vector containing the number of Floquet
  %           points in the infinite directions.
  %
  %         * opts, structure with fields
  %           -> solBasis, if set to true, allows to compute the solution
  %              for any of the basis functions;
  %           -> omega, the frequency (needed to compute the flux);
  %           -> suffix for post-processing.
  %
  % OUTPUTS: * U, a structure containing the solution in each cell of periodicity.

  %% % ************************* %
  %  % Preliminary verifications %
  %  % ************************* %
  pgbars = findall(0,'type','figure','tag','TMWWaitbar');
  delete(pgbars);

  infiniteDirections = unique(infiniteDirections(:).');
  phi = BCstruct.phi;
  opts.verbose = 0;

  if (abs(orientation) ~= 1)
    % Produce an error if an inappropriate value is given for the
    % sign of the outward normal
    error('%d a ete donne comme orientation du demi-espace au lieu de -1 ou de 1.', orientation);
  end

  if min((1:mesh.dimension) ~= semiInfiniteDirection)
    % The component along which the domain is semi-infinite should be between
    % 1 and the mesh dimension
    error('Le domaine ne peut pas être semi-infini dans la direction %d.', semiInfiniteDirection);
  end

  for idI = 1:length(infiniteDirections)
    if min((1:mesh.dimension) ~= infiniteDirections(idI))
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
  N = mesh.numPoints;
  Ni = length(infiniteDirections);

  AAelem = cell(Ni +1, Ni + 1);
  for idU = 1:(Ni + 1)
    for idV = 1:(Ni + 1)
      AAelem{idV, idU} = sparse(N, N);
    end
  end

  for idT = 1:length(volBilinearIntg.alpha_u)
    fun = volBilinearIntg.fun{idT};

    Au = cell(Ni + 1, 1);
    Av = cell(Ni + 1, 1);
    Au{1} = volBilinearIntg.alpha_u{idT};
    Av{1} = volBilinearIntg.alpha_v{idT};
    for idB = 1:Ni
      Au{idB + 1} = Au{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
      Av{idB + 1} = Av{1}(:, infiniteDirections(idB) + 1) * [1 0 0 0];
    end

    for idU = 1:(Ni + 1)
      for idV = 1:(Ni + 1)
        AAelem{idV, idU} = AAelem{idV, idU} + FEPack.pdes.Form.global_matrix(mesh.domain('volumic'), Au{idU}, Av{idV}, fun);
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
    AA = sparse(N, N);
    for idU = 1:(Ni + 1)
      for idV = 1:(Ni + 1)
        AA = AA + FBvar(idV)' * AAelem{idV, idU} * FBvar(idU);
      end
    end

    % The Floquet-Bloch transform of the boundary data along the relevant directions
    BCstruct.phi = @(x) FEPack.tools.BlochTransform(x, K(idFB, :), phi, infiniteDirections);

    % Compute the solution
    TFBU{idFB} = PeriodicHalfGuideBVP(mesh, orientation, semiInfiniteDirection, AA, BCstruct, numCellsSemiInfinite, opts);

  end

  delete(pgbar);

  %% % ***************************************** %
  %  % Apply the inverse Floquet-Bloch transform %
  %  % ***************************************** %
  numCells = [1 1 1];
  numCells(infiniteDirections) = 2*numCellsInfinite;
  numCells(semiInfiniteDirection) = numCellsSemiInfinite;
  Nu = prod(numCells);

  [I1, I2, I3] = ind2sub(numCells, 1:Nu);
  pointsIds = [I1; I2; I3]; % 3-by-Nu
  tau = pointsIds(infiniteDirections, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
  W = prod(2*pi ./ (numFloquetPoints - 1)); % 1-by-1
  U = zeros(N, Nu);

  for idFB = 1:numFBpoints
    % The integral that defines the inverse Floquet-Bloch transform is computed
    % using a rectangular rule.
    exp_k_dot_x = exp(1i * mesh.points(:, infiniteDirections) * K(idFB, :).'); % N-by-1
    exp_k_dot_tau = exp(1i * K(idFB, :) * tau); % 1-by-Nu
    U_TFB = TFBU{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

    U = U + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
  end

  U = U / sqrt(2*pi);

end
