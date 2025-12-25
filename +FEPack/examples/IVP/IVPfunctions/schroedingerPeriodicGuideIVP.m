function out = schroedingerPeriodicGuideIVP(...
  delta_t, final_t,...
  infiniteDirection,...
  AApos, BBpos, mshPos, BCstructPos, numCellsPos,...
  AAneg, BBneg, mshNeg, BCstructNeg, numCellsNeg,...
  AAint, BBint, vecInitial, mshInt, spBint_pos, spBint_neg,...
  rhsTheta)

  % schroedingerPeriodicGuideIVP

  %  ******************************* %
  %% Initialize variables and output %
  %  ******************************* %
  numTsteps = ceil(final_t / delta_t);

  Npos  = mshPos.numPoints;
  NbPos = BCstructPos.spB0.numBasis;

  tmp.pos.prop = cell(numTsteps, 1);
  tmp.pos.DtN  = cell(numTsteps, 1);
  out.pos.sol  = zeros(Npos, numCellsPos, numTsteps);
  auxSolPos    = zeros(Npos, NbPos, numCellsPos, numTsteps);

  Nneg  = mshNeg.numPoints;
  NbNeg = BCstructNeg.spB0.numBasis;

  tmp.neg.prop = cell(numTsteps, 1);
  tmp.neg.DtN  = cell(numTsteps, 1);
  out.neg.sol  = zeros(Nneg, numCellsNeg, numTsteps);
  auxSolNeg    = zeros(Nneg, NbNeg, numCellsNeg, numTsteps);

  out.int.sol = zeros(mshInt.numPoints, numTsteps);

  %  *************** %
  %% First iteration %
  %  *************** %
  %% Positive side: auxiliary half-guide solution
  [~, E1pos, ~, T10pos, ~, T11pos,...
  tmp.pos.prop{1}, tmp.pos.DtN{1}, auxSolPos(:, :, :, 1)] =...
  halfGuideDirichlet(infiniteDirection, mshPos, BCstructPos, AApos, numCellsPos);

  T11_p_DtN_inv_pos = (T11pos + tmp.pos.DtN{1}) \ eye(size(T11pos));
  % theta = 0.5;
  % co = @(n) sqrt((2*pi*n).^2 - 1i / (delta_t*theta));

  %% Negative side: auxiliary half-guide solution
  [~, E1neg, ~, T10neg, ~, T11neg,...
  tmp.neg.prop{1}, tmp.neg.DtN{1}, auxSolNeg(:, :, :, 1)] =...
  halfGuideDirichlet(infiniteDirection, mshNeg, BCstructNeg, AAneg, numCellsNeg);

  T11_p_DtN_inv_neg = (T11neg + tmp.neg.DtN{1}) \ eye(size(T11neg));

  %% Problem in interior domain
  % Find bounded directions and impose periodic condition on them
  boundedDirections = find((1:mshInt.dimension) ~= infiniteDirection); % Bounded directions

  u = FEPack.pdes.PDEObject;
  ecsInt = FEPack.pdes.EssentialConditions;
  for idI = 1:length(boundedDirections)
    ecsInt = ecsInt & assignEcs((u|mshInt.domains{2*boundedDirections(idI)-1}) -...
                                (u|mshInt.domains{2*boundedDirections(idI)}), 0.0);
  end

  % Domains
  Sigma_neg = mshInt.domains{2*infiniteDirection};
  Sigma_pos = mshInt.domains{2*infiniteDirection-1};

  % Coefficients associated to boundary condition
  SSpos = FEPack.pdes.Form.intg_TU_V(Sigma_pos, tmp.pos.DtN{1}, BCstructPos.spB0, 'projection');
  SSneg = FEPack.pdes.Form.intg_TU_V(Sigma_neg, tmp.neg.DtN{1}, BCstructNeg.spB0, 'projection');

  AAint = AAint + SSpos + SSneg;

  % Right-hand side
  LLint = BBint * vecInitial - rhsTheta(1);

  % Solve interior problem
  out.int.sol(:, 1) = CellBVP(mshInt, AAint, LLint, ecsInt);

  %% Construct the solution in the whole domain
  % Positive side
  phi_pos = spBint_pos.FE_to_spectral * out.int.sol(Sigma_pos.IdPoints, 1);
  for idCell = 1:numCellsPos
    out.pos.sol(:, idCell, 1) = auxSolPos(:, :, idCell, 1) * phi_pos;
  end

  % Negative side
  phi_neg = spBint_neg.FE_to_spectral * out.int.sol(Sigma_neg.IdPoints, 1);
  for idCell = 1:numCellsNeg
    out.neg.sol(:, idCell, 1) = auxSolNeg(:, :, idCell, 1) * phi_neg;
  end

  %  *************** %
  %% Next iterations %
  %  *************** %
  % Auxiliary solutions and trace operators
  [auxEpos, auxF0pos, auxF1pos] = auxiliary_cell_solution(infiniteDirection, mshPos, BCstructPos, AApos, BBpos);

  [auxEneg, auxF0neg, auxF1neg] = auxiliary_cell_solution(infiniteDirection, mshNeg, BCstructNeg, AAneg, BBneg);

  fixed_point.numIter = 1e3;
  fixed_point.tol     = 1e-4;

  for idT = 2:numTsteps

    %% Positive side: current propagation operator
    % Fixed-point iteration to solve the generalized Sylvester equation
    SylvRhsPos = auxF1pos * auxSolPos(:, :, 1, idT-1) +...
                auxF0pos * auxSolPos(:, :, 1, idT-1) * tmp.pos.prop{1};
    for idL = 2:idT-1
      SylvRhsPos = SylvRhsPos + tmp.pos.DtN{idL} * tmp.pos.prop{idT+1-idL};
    end

    X = tmp.pos.prop{idT-1};
    idI = 0;
    while (true)
      Xnew = -T11_p_DtN_inv_pos * (SylvRhsPos + T10pos * X * tmp.pos.prop{1});
      ecart = norm(Xnew + T11_p_DtN_inv_pos * (SylvRhsPos + T10pos * Xnew * tmp.pos.prop{1}), 'fro') / norm(Xnew, 'fro');
      
      if (ecart < fixed_point.tol)
        fprintf('%d: Fixed-point iteration tolerance reached in %d iteration(s).\n', idT, idI+1);
        break;
      elseif (idI >= fixed_point.numIter)
        fprintf('%d: Fixed-point maximum iteration reached with error %0.5e instead of %0.5e.\n', idT, ecart, fixed_point.tol);
        break;
      else
        idI = idI + 1;
        X = Xnew;
      end
    end

    % DtN operator
    tmp.pos.prop{idT} = Xnew;
    tmp.pos.DtN{idT}  = auxF0pos * auxSolPos(:, :, 1, idT-1) + T10pos * tmp.pos.prop{idT};

    %% Positive side: cell by cell construction
    auxSolPos(:, :, 1, idT) = auxEpos * auxSolPos(:, :, 1, idT-1) + E1pos * tmp.pos.prop{idT};

    for idCell = 2:numCellsPos
      for idL = 1:idT
        auxSolPos(:, :, idCell, idT) = auxSolPos(:, :, idCell, idT) +...
          auxSolPos(:, :, idCell-1, idL) * tmp.pos.prop{idT-idL+1};
      end
    end

    %% Negative side: current propagation operator
    % Fixed-point iteration to solve the generalized Sylvester equation
    SylvRhsNeg = auxF1neg * auxSolNeg(:, :, 1, idT-1) +...
                auxF0neg * auxSolNeg(:, :, 1, idT-1) * tmp.neg.prop{1};
    for idL = 2:idT-1
      SylvRhsNeg = SylvRhsNeg + tmp.neg.DtN{idL} * tmp.neg.prop{idT+1-idL};
    end

    X = tmp.neg.prop{idT-1};
    idI = 0;
    while (true)
      Xnew = -T11_p_DtN_inv_neg * (SylvRhsNeg + T10neg * X * tmp.neg.prop{1});
      ecart = norm(Xnew + T11_p_DtN_inv_neg * (SylvRhsNeg + T10neg * Xnew * tmp.neg.prop{1}), 'fro') / norm(Xnew, 'fro');

      if (ecart < fixed_point.tol)
        fprintf('%d: Fixed-point iteration tolerance reached in %d iteration(s).\n', idT, idI+1);
        break;
      elseif  (idI >= fixed_point.numIter)
        fprintf('%d: Fixed-point maximum iteration reached with error %0.5e instead of %0.5e.\n', idT, ecart, fixed_point.tol)
        break;
      else
        idI = idI + 1;
        X = Xnew;
      end
    end

    % DtN operator
    tmp.neg.prop{idT} = Xnew;
    tmp.neg.DtN{idT}  = auxF0neg * auxSolNeg(:, :, 1, idT-1) + T10neg * tmp.neg.prop{idT};

    %% Negative side: cell by cell construction
    auxSolNeg(:, :, 1, idT) = auxEneg * auxSolNeg(:, :, 1, idT-1) + E1neg * tmp.neg.prop{idT};

    for idCell = 2:numCellsNeg
      for idL = 1:idT
        auxSolNeg(:, :, idCell, idT) = auxSolNeg(:, :, idCell, idT) +...
          auxSolNeg(:, :, idCell-1, idL) * tmp.neg.prop{idT-idL+1};
      end
    end
    
    %% Interior problem
    LLint = BBint * out.int.sol(:, idT-1) - rhsTheta(idT);

    % Surface contributions
    for idL = 2:idT
      phi_pos = zeros(mshInt.numPoints, 1);
      phi_pos(Sigma_pos.IdPoints) = out.int.sol(Sigma_pos.IdPoints, idT-idL+1);
      SSpos = FEPack.pdes.Form.intg_TU_V(Sigma_pos, tmp.pos.DtN{idL}, BCstructPos.spB0, 'projection');

      phi_neg = zeros(mshInt.numPoints, 1);
      phi_neg(Sigma_neg.IdPoints) = out.int.sol(Sigma_neg.IdPoints, idT-idL+1);
      SSneg = FEPack.pdes.Form.intg_TU_V(Sigma_neg, tmp.neg.DtN{idL}, BCstructNeg.spB0, 'projection');

      LLint = LLint - SSpos * phi_pos - SSneg * phi_neg;
    end

    % Solve interior problem
    out.int.sol(:, idT) = CellBVP(mshInt, AAint, LLint, ecsInt);

    %% Construct entire solution
    for idCell = 1:numCellsPos
      for idL = 1:idT
        phi_pos = spBint_pos.FE_to_spectral * out.int.sol(Sigma_pos.IdPoints, idT-idL+1);

        out.pos.sol(:, idCell, idT) = out.pos.sol(:, idCell, idT) +...
          auxSolPos(:, :, idCell, idL) * phi_pos;
      end
    end
    
    for idCell = 1:numCellsNeg
      for idL = 1:idT
        phi_neg = spBint_neg.FE_to_spectral * out.int.sol(Sigma_neg.IdPoints, idT-idL+1);

        out.neg.sol(:, idCell, idT) = out.neg.sol(:, idCell, idT) +...
          auxSolNeg(:, :, idCell, idL) * phi_neg;
      end
    end
  end

end

function [E0, E1, T00, T10, T01, T11, R, DtN, U] = halfGuideDirichlet(infiniteDirection, msh, BCstruct, AA, numCells)

  %% Initialization
  N = msh.numPoints;
  
  % Find bounded directions and impose periodic condition on them
  boundedDirections = find((1:msh.dimension) ~= infiniteDirection); % Bounded directions

  u = FEPack.pdes.PDEObject;
  ecs = FEPack.pdes.EssentialConditions;
  for idI = 1:length(boundedDirections)
    ecs = ecs & assignEcs((u|msh.domains{2*boundedDirections(idI)-1}) - (u|msh.domains{2*boundedDirections(idI)}), 0.0);
  end

  % Domains
  Sigma0 = msh.domains{2*infiniteDirection};
  Sigma1 = msh.domains{2*infiniteDirection-1};

  % Boundary basis functions
  spB0 = BCstruct.spB0;
  spB1 = BCstruct.spB1;
  Nb = spB0.numBasis; % Number of basis functions

  %% Solve the local cell problems
  % Periodic cell problems whose boundary conditions are similar to the one of
  % the half-guide problem

  % Surfacic right-hand sides
  B0 = sparse(N, Nb);
  if (spB0.is_interpolated)
    B0(Sigma0.IdPoints, :) = spB0.phis;
  else
    B0(Sigma0.IdPoints, :) = spB0.phis(msh.points(Sigma0.IdPoints, :), 1:Nb);
  end

  B1 = sparse(N, Nb);
  if (spB1.is_interpolated)
    B1(Sigma1.IdPoints, :) = spB1.phis;
  else
    B1(Sigma1.IdPoints, :) = spB1.phis(msh.points(Sigma1.IdPoints, :), 1:Nb);
  end

  BCstruct.representation = 'projection';

  % Homogeneous Dirichlet condition
  ecs = ecs & assignEcs(u|Sigma0, 0.0) & assignEcs(u|Sigma1, 0.0);
  ecs.applyEcs;

  % Add the surfacic contributions
  ecs.b = [B0, B1];

  % Solve the local cell problems
  AA0 =  ecs.P * AA * ecs.P';
  LL0 = -ecs.P * AA * ecs.b;

  tic;
  UU0 = AA0 \ LL0;
  toc;

  Ecell = ecs.b + ecs.P' * UU0;

  % Deduce the local cell solutions.
  E0 = Ecell(:, 1:Nb);
  E1 = Ecell(:, 1+Nb:end);

  % Traces and normal traces of the local cell solutions
  % Ekl is the trace of Ek on Sigmal
  E00 = speye(Nb);
  E10 = sparse(Nb, Nb);
  E01 = sparse(Nb, Nb);
  E11 = speye(Nb);

  % Fkl is the normal trace of Ek on Sigmal. It is evaluated weakly
  AAE0 = AA * E0;
  AAE1 = AA * E1;

  MME0 = spB0.invmassmat * B0';
  MME1 = spB0.invmassmat * B1';

  % NE COMPTE PAS ICI CAR ON A MULTIPLIE PAR
  % SPB0.INVMASSMAT
  % La multiplication par la variable + ou - 1       %
  % orientation est nécessaire, contrairement à la   % 
  % théorie, car dans le cas négatif, le produit     %
  %       (Ej' * AA * El)                            % 
  % est égal à une double intégrale avec la variable %
  % semi-infinie qui ***part de 0 à -perNeg***, à    %
  % cause de la façon dont mon maillage est défini.  %
  % ************************************************ %
  T00 = MME0 * AAE0;
  T10 = MME0 * AAE1;
  T01 = MME1 * AAE0;
  T11 = MME1 * AAE1;

  %% Solve the Riccati equation
  flux = @(V) error('Tu ne devrais pas arriver ici.');
  riccatiOpts.tol = 1.0e-2;
  riccatiOpts.suffix = '';
  R = propagationOperators([E01, E11;  T01,  T11], ...
                           [E00, E10; -T00, -T10], flux, riccatiOpts);

  %% Compute the solution cell by cell
  U = zeros(N, Nb, numCells);
  Rn = eye(Nb);

  for idCell = 0:numCells-1
    % Compute the solution in the current cell
    U(:, :, idCell + 1) = E0 * Rn + E1 * R * Rn;

    % Update
    Rn = R * Rn;
  end
  
  %% The DtN operator
  DtN = T00 + T10 * R;
  
end


function [auxE, auxF0, auxF1] = auxiliary_cell_solution(infiniteDirection, msh, BCstruct, AA, BB)

  % Find bounded directions and impose periodic condition on them
  boundedDirections = find((1:msh.dimension) ~= infiniteDirection); % Bounded directions

  u = FEPack.pdes.PDEObject;
  ecs = FEPack.pdes.EssentialConditions;
  for idI = 1:length(boundedDirections)
    ecs = ecs & assignEcs((u|msh.domains{2*boundedDirections(idI)-1}) - (u|msh.domains{2*boundedDirections(idI)}), 0.0);
  end

  % Domains
  Sigma0 = msh.domains{2*infiniteDirection};
  Sigma1 = msh.domains{2*infiniteDirection-1};

  % Boundary basis functions
  spB0 = BCstruct.spB0;
  spB1 = BCstruct.spB1;

  % Homogeneous Dirichlet condition
  ecs = ecs & assignEcs(u|Sigma0, 0.0) & assignEcs(u|Sigma1, 0.0);
  ecs.applyEcs;

  %% Solve the local cell problems
  AA0 = ecs.P * AA * ecs.P';
  BB0 = ecs.P * BB;

  tic;
  auxE = ecs.P' * (AA0 \ BB0);
  toc;
  
  %% Surface operators
  B0 = sparse(msh.numPoints, spB0.numBasis);
  if (spB0.is_interpolated)
    B0(Sigma0.IdPoints, :) = spB0.phis;
  else
    B0(Sigma0.IdPoints, :) = spB0.phis(msh.points(Sigma0.IdPoints, :), 1:spB0.numBasis);
  end

  B1 = sparse(msh.numPoints, spB1.numBasis);
  if (spB1.is_interpolated)
    B1(Sigma1.IdPoints, :) = spB1.phis;
  else
    B1(Sigma1.IdPoints, :) = spB1.phis(msh.points(Sigma1.IdPoints, :), 1:spB1.numBasis);
  end

  % NE COMPTE PAS ICI CAR ON A MULTIPLIE PAR
  % SPB0.INVMASSMAT
  % Normal traces of the auxiliary solution
  % Here again, it is very important to multiply by
  % orientation, because of how my mesh is defined on
  % the negative domain.
  auxF0 = spB0.invmassmat * B0' * (AA * auxE - BB);
  auxF1 = spB1.invmassmat * B1' * (AA * auxE - BB);

end