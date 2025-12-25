function U = PeriodicGuideBVP(infiniteDirection,...
  volBilinearIntg_pos, mesh_pos, BCstruct_pos, numCells_pos,...
  volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCells_neg,...
  volBilinearIntg_int, volLinearIntg, mesh_int, spBint_pos, spBint_neg, opts)

  % PeriodicGuideBVP
  if (nargin < 14)
    opts = [];
  end
  
  % Preliminary set ups
  opts.solBasis = true;

  % Solve the problem in the positive half-guide
  % ////////////////////////////////////////////
  [Upos, ~, ~, ~, ~, BCstruct_pos, Lambda_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, infiniteDirection, volBilinearIntg_pos, BCstruct_pos, numCells_pos, opts);

  % Solve the problem in the negative half-guide
  % ////////////////////////////////////////////
  [Uneg, ~, ~, ~, ~, BCstruct_neg, Lambda_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, infiniteDirection, volBilinearIntg_neg, BCstruct_neg, numCells_neg, opts);

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % N = 32;
  % gamma = sqrt(4*pi*pi*(-N/4:N/4).^2 - opts.omega^2);
  % Lambda_pos = diag(gamma);
  % Lambda_neg = diag(gamma);
  % Nbpos = BCstruct_pos.spB0.numBasis;
  % Nbneg = BCstruct_neg.spB0.numBasis;
  %
  % for i = 1:numCells_pos
  %   Upos{i} = BCstruct_pos.spB0.phis(mesh_pos.points, 1:Nbpos) .* exp(-(mesh_pos.points(:, 1) + (i-1))*gamma);
  % end
  % for i = 1:numCells_neg
  %   Uneg{i} = BCstruct_neg.spB0.phis(mesh_neg.points, 1:Nbneg) .* exp((mesh_neg.points(:, 1) - (i-1))*gamma);
  % end
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Solve the problem in the interior domain
  % ////////////////////////////////////////
  % volumic FE matrix and right-hand side
  if isa(volBilinearIntg_int, 'FEPack.pdes.Form')
    % Compute the FE matrix if not done already
    AAint = FEPack.pdes.Form.intg(mesh_int.domain('volumic'), volBilinearIntg_int);
  else
    AAint = volBilinearIntg_int;
  end

  if isa(volLinearIntg, 'FEPack.pdes.Form')
    % Compute the FE right-hand side if not done already
    LLint = FEPack.pdes.Form.intg(mesh_int.domain('volumic'), volLinearIntg);
  else
    LLint = volLinearIntg;
  end

  % Find bounded directions and impose periodic condition on them
  boundedDirections = find((1:mesh_int.dimension) ~= infiniteDirection); % Bounded directions

  u = FEPack.pdes.PDEObject;
  ecs = FEPack.pdes.EssentialConditions;
  for idI = 1:length(boundedDirections)
    ecs = ecs & assignEcs((u|mesh_int.domains{2*boundedDirections(idI)-1}) -...
                          (u|mesh_int.domains{2*boundedDirections(idI)}), 0.0);
  end

  % Domains
  Sigma_neg = mesh_int.domains{2*infiniteDirection};
  Sigma_pos = mesh_int.domains{2*infiniteDirection-1};

  % Coefficients associated to boundary condition
  Nbpos = size(Lambda_pos, 1);
  Nbneg = size(Lambda_neg, 1);
  BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * Lambda_pos) \ (Lambda_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
  BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * Lambda_neg) \ (Lambda_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

  dirichletBC = (norm(BCu_int_pos, 'fro') < eps) || (norm(BCu_int_neg, 'fro') < eps);

  if (dirichletBC)
    error(['Le problème de bord semble contenir une condition ',...
           'de Dirichlet. Revoir les problèmes de demi-guide considérés.'])
  end

  SSpos = FEPack.pdes.Form.intg_TU_V(Sigma_pos, BCu_int_pos, BCstruct_pos.spB0, 'projection');
  SSneg = FEPack.pdes.Form.intg_TU_V(Sigma_neg, BCu_int_neg, BCstruct_neg.spB0, 'projection');

  AAint = AAint - SSpos - SSneg;

  % Solve interior problem
  U.interior = CellBVP(mesh_int, AAint, LLint, ecs);

  % Construct the solution in the whole domain
  % //////////////////////////////////////////
  % Positive side
  phi_pos = (BCstruct_pos.BCu * eye(Nbpos) - BCstruct_pos.BCdu * BCu_int_pos) * spBint_pos.FE_to_spectral * U.interior(Sigma_pos.IdPoints);
  U.positive = zeros(mesh_pos.numPoints, numCells_pos);

  for idCell = 1:numCells_pos
    U.positive(:, idCell) = Upos{idCell} * phi_pos;
  end

  % Negative side
  phi_neg = (BCstruct_neg.BCu * eye(Nbneg) - BCstruct_neg.BCdu * BCu_int_neg) * spBint_neg.FE_to_spectral * U.interior(Sigma_neg.IdPoints);
  U.negative = zeros(mesh_neg.numPoints, numCells_neg);

  for idCell = 1:numCells_neg
    U.negative(:, idCell) = Uneg{idCell} * phi_neg;
  end
end
