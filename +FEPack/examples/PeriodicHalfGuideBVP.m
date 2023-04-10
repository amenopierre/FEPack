 function [U, newBCstruct, Lambda] = PeriodicHalfGuideBVP(mesh, orientation, infiniteDirection, volBilinearIntg, BCstruct, numCells, opts)
  % PeriodicHalfGuideBVP(mesh, orientation, infiniteDirection, volBilinearIntg, BCstruct, numCells, opts) solves
  %
  %       Find u in V such that tau(u) = g, and
  %                A(u, v) + S(u, v) = Sf(v) for any v in V such that tau(v) = 0,
  %
  % where tau is a 0-order linear operator.
  %
  % INPUTS: * mesh, Mesh object.
  %
  %         * orientation, a number in {-1, 1} that represent whether the guide
  %           is infinite along a positive or a negative axis.
  %
  %         * infiniteDirection, direction along which the guide is infinite.
  %           infiniteDirection is between 1 and the mesh dimension.
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
  %         * numCells, Number of cells on which the solution is computed.
  %
  %         * opts, structure with fields
  %           -> solBasis, if set to true, allows to compute the solution
  %              for any of the basis functions;
  %           -> omega, the frequency (needed to compute the flux);
  %           -> verbose, if set to 0, prevents the code from printing any
  %              superfluous message;
  %           -> suffix for post-processing.
  %
  % OUTPUTS: * U, a structure containing the solution in each cell of periodicity.
  %              -> if opts.solBasis = true, then U is a numCells-sized cell
  %                 array. Each component of U is a N-by-Nb matrix, where Nb is
  %                 the number of basis functions.
  %              -> if opts.solBasis = false, then U is a N-by-numCells matrix.
  %
  %          * newBCstruct, BCstruct with some modified fields.
  %
  %          * Lambda, linear combination of the trace and the normal trace
  %            of U on the boundary {x_i = 0} (if the guide is infinite along x_i).

  %% % ************************* %
  %  % Preliminary verifications %
  %  % ************************* %
  if (abs(orientation) ~= 1)
    % Produce an error if an inappropriate value is given for the
    % sign of the outward normal
    error('%d a ete donne comme orientation du demi-guide au lieu de -1 ou de 1.', orientation);
  end

  if min((1:mesh.dimension) ~= infiniteDirection)
    % The component along which the guide is infinite should be between
    % 1 and the mesh dimension
    error('Le demi-guide ne peut pas être infini dans la direction %d.', infiniteDirection);
  end

  if ~isfield(opts, 'solBasis')
    opts.solBasis = false;
  end

  if ~isfield(opts, 'verbose')
    opts.verbose = 1;
  end

  %% % ************** %
  %  % Initialization %
  %  % ************** %
  N = mesh.numPoints;
  if isa(volBilinearIntg, 'FEPack.pdes.Form')
    % Compute the FE matrix if not done already
    AA = FEPack.pdes.Form.intg(mesh.domain('volumic'), volBilinearIntg);
  else
    AA = volBilinearIntg;
  end

  % Find bounded directions and impose periodic condition on them
  boundedDirections = find((1:mesh.dimension) ~= infiniteDirection); % Bounded directions

  u = FEPack.pdes.PDEObject;
  v = dual(u);
  ecs = FEPack.pdes.EssentialConditions;
  for idI = 1:length(boundedDirections)
    ecs = ecs & assignEcs((u|mesh.domains{2*boundedDirections(idI)-1}) - (u|mesh.domains{2*boundedDirections(idI)}), 0.0);
  end

  % Domains
  Sigma0 = mesh.domains{2*infiniteDirection};
  Sigma1 = mesh.domains{2*infiniteDirection-1};

  % Boundary basis functions
  spB0 = BCstruct.spB0;
  spB1 = BCstruct.spB1;
  Nb = spB0.numBasis; % Number of basis functions

  %% % ***************************** %
  %  % Solve the local cell problems %
  %  % ***************************** %
  if (opts.verbose)
    fprintf('1. Resolution des problemes de cellule locaux\n');
  end
  % Periodic cell problems whose boundary conditions are similar to the one of
  % the half-guide problem

  % Surfacic right-hand sides
  B0 = sparse(N, Nb);
  if (spB0.is_interpolated)
    B0(Sigma0.IdPoints, :) = spB0.phis;
  else
    B0(Sigma0.IdPoints, :) = spB0.phis(mesh.points(Sigma0.IdPoints, :), 1:Nb);
  end

  B1 = sparse(N, Nb);
  if (spB1.is_interpolated)
    B1(Sigma1.IdPoints, :) = spB1.phis;
  else
    B1(Sigma1.IdPoints, :) = spB1.phis(mesh.points(Sigma1.IdPoints, :), 1:Nb);
  end

  if (abs(BCstruct.BCdu) < eps)

    % Dirichlet boundary conditions %
    % ***************************** %
    if (opts.verbose)
      fprintf('\tConditions de Dirichlet\n');
    end

    if (isfield(BCstruct, 'BCu') && ...
      (~isa(BCstruct.BCu, 'double') || (BCstruct.BCu ~= 1.0)))
      warning('on');
      warning(['Pour des conditions de Dirichlet, le coefficient devant u ',...
               'dans la condition aux bords est automatiquement pris égal à 1.']);
    end
    BCstruct.BCu = 1.0;
    BCstruct.representation = 'projection';

    % Homogeneous Dirichlet condition
    ecs = ecs & assignEcs(u|Sigma0, 0.0) & assignEcs(u|Sigma1, 0.0);
    ecs.applyEcs;

    % Add the surfacic contributions
    ecs.b = [B0, B1];

    % Solve the local cell problems
    Ecell = CellBVP(mesh, AA, 0.0, ecs);

    % Deduce the local cell solutions.
    E0 = Ecell(:, 1:Nb);
    E1 = Ecell; E1(:, 1:Nb) = [];

    if (opts.verbose)
      fprintf('2. Traces et traces normales\n');
    end

    % Traces and normal traces of the local cell solutions
    % Ekl is the trace of Ek on Sigmal
    E00 = speye(Nb);
    E10 = sparse(Nb, Nb);
    E01 = sparse(Nb, Nb);
    E11 = speye(Nb);

    % Fkl is the normal trace of Ek on Sigmal. It is evaluated weakly
    F00 = spB0.invmassmat * (E0' * AA * E0);
    F10 = spB0.invmassmat * (E0' * AA * E1);
    F01 = spB0.invmassmat * (E1' * AA * E0);
    F11 = spB0.invmassmat * (E1' * AA * E1);

  else

    % Robin (Neumann) boundary conditions %
    % *********************************** %
    if (opts.verbose)
      fprintf('\tConditions de Robin-Neumann\n');
    end

    % Matrix associated to the multiplication by a function
    if isa(BCstruct.BCu, 'function_handle')
      BCstruct.BCu = spB0.intg_U_V(Sigma0, BCstruct.BCu);
      BCstruct.representation = 'weak evaluation';
    end

    % Surfacic contributions
    if (~isfield(BCstruct, 'representation') || (length(BCstruct.BCu) == 1))
      BCstruct.representation = 'projection';
    end

    SS0 = FEPack.pdes.Form.intg_TU_V(Sigma0, BCstruct.BCu, BCstruct.representation);
    SS1 = FEPack.pdes.Form.intg_TU_V(Sigma1, BCstruct.BCu, BCstruct.representation);
    AA = AA + SS0 + SS1;

    % Surfacic right-hand sides
    BB = [FEPack.pdes.Form.intg(Sigma0, id(u)*id(v)) * B0,...
          FEPack.pdes.Form.intg(Sigma1, id(u)*id(v)) * B1];

    % Solve the local cell problems
    Ecell = CellBVP(mesh, AA, BB, ecs);

    % Deduce the local cell solutions.
    E0 = Ecell(:, 1:Nb);
    E1 = Ecell; E1(:, 1:Nb) = [];

    if (opts.verbose)
      fprintf('2. Traces et traces normales\n');
    end

    % Ekl is the trace of Ek on Sigmal
    E00 = spB0.FE_to_spectral * E0(Sigma0.IdPoints, :);
    E10 = spB1.FE_to_spectral * E1(Sigma0.IdPoints, :);
    E01 = spB0.FE_to_spectral * E0(Sigma1.IdPoints, :);
    E11 = spB1.FE_to_spectral * E1(Sigma1.IdPoints, :);

    % Set BCu to projection matrix if not done already
    if strcmpi(BCstruct.representation, 'weak evaluation')
      BCstruct.BCu = BCstruct.spB0.invmassmat * BCstruct.BCu;
      BCstruct.representation = 'projection';
    end

    % Fkl is the normal trace of Ek on Sigmal.
    invBCdu = 1.0 / BCstruct.BCdu;
    F00 = -invBCdu * (BCstruct.BCu * E00 - speye(Nb));
    F10 = -invBCdu *  BCstruct.BCu * E10;
    F01 = -invBCdu *  BCstruct.BCu * E01;
    F11 = -invBCdu * (BCstruct.BCu * E11 - speye(Nb));

  end

  %% % ************************** %
  %  % Solve the Riccati equation %
  %  % ************************** %
  if (opts.verbose)
    fprintf('3. Resolution du système de Riccati\n');
  end

  % Solve the linearized eigenvalue problem associated to the Riccati equation
  flux = @(V) modesFlux(V, E00, E10, F00, F10, orientation, spB0, opts.omega);
  riccatiOpts.tol = 1.0e-2;
  if isfield(opts, 'suffix')
    riccatiOpts.suffix = opts.suffix;
  else
    riccatiOpts.suffix = '';
  end
  [R, D] = propagationOperators([E01, E11;  orientation*F01,  orientation*F11], ...
                                [E00, E10; -orientation*F00, -orientation*F10], flux, riccatiOpts);

  %% % ********************************* %
  %  % Compute the solution cell by cell %
  %  % ********************************* %
  if (opts.verbose)
    fprintf('4. Construction de la solution\n');
  end

  if (opts.solBasis)

    % Compute the solution for any basis function
    U = cell(numCells, 1);
    R0 = eye(Nb);
    R1 = D;

    for idCell = 0:numCells-1
      % Compute the solution in the current cell
      U{idCell + 1} = E0 * R0 + E1 * R1;

      % Update the coefficients
      R0 = R * R0;
      R1 = D * R0;
    end

  else

    % Compute the solution for one boundary data
    U = zeros(N, numCells);
    R0Phi = spB0.FE_to_spectral * BCstruct.phi(mesh.points(Sigma0.IdPoints, :));
    R1Phi = D * R0Phi;

    for idCell = 0:numCells-1
      % Compute the solution in the current cell
      U(:, idCell + 1) = E0 * R0Phi + E1 * R1Phi;

      % Update the coefficients
      R0Phi = R * R0Phi;
      R1Phi = D * R0Phi;
    end

  end

  %% % *********************************************** %
  %  % The transmission coefficient and the derivative %
  %  % *********************************************** %
  if (opts.verbose)
    fprintf('5. Calcul de l''operateur de transmission global\n');
  end

  % Compute the trace and the normal trace of the half-guide solution
  U0 = E00 + E10 * D;
  dU0 = F00 + F10 * D;

  Lambda = -BCstruct.BCu' * dU0 + BCstruct.BCdu' * U0;
  newBCstruct = BCstruct;
end

%% % ******** %
%  % Appendix %
%  % ******** %
% The flux function
function flux = modesFlux(V, E00, E10, F00, F10, orientation, spB, omega)
  
  Nb = size(E00, 1);
  Psi  = E00 * V(1:Nb, :) + E10 * V(Nb+1:end, :);
  dPsi = F00 * V(1:Nb, :) + F10 * V(Nb+1:end, :);

  flux = -orientation * imag(diag(Psi' * spB.massmat * dPsi) / omega);

end
