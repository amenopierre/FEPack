function [dispmesh, val] = edge_states_irrational_3D(...
  BBkpar, numNodesK,...
  BBlambda, numNodesLambda,...
  op_pos, numNodesXpos, numNodesYpos, numNodesZpos, bc_basis_functions_pos,...
  op_neg, numNodesXneg, numNodesYneg, numNodesZneg, bc_basis_functions_neg,...
  op_int, BBint, numNodesXint, numNodesYint, numNodesZint,...
  parallel_use, parallel_numcores)
  
  %% [mesh, val] = EDGE_STATES_RATIONAL()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Meshes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% Mesh associated to parallel quasi-momentum / spectral parameter
  dispmesh = FEPack.meshes.MeshRectangle(0, BBkpar, BBlambda, numNodesK, numNodesLambda);

  %% ========= (+) half-guide
  % Mesh
  mesh_pos = FEPack.meshes.MeshCuboid(1, [0  1], [0 1], [0 1], numNodesXpos, numNodesYpos, numNodesZpos);
  cell_pos = mesh_pos.domain('volumic');

  % Boundary conditions
  if strcmpi(bc_basis_functions_pos, 'Lagrange')
    init_BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domain('xmin'));
    init_BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domain('xmax'));
  else
    FourierIds = [0 0 0]; FourierIds(1) = floor(numNodesXpos/4);
    init_BCstruct_pos.spB0 = FEPack.spaces.FourierBasis(mesh_pos.domain('xmin'), FourierIds);
    init_BCstruct_pos.spB1 = FEPack.spaces.FourierBasis(mesh_pos.domain('xmax'), FourierIds);
  end

  %% ========= (-) half-guide
  % Mesh
  mesh_neg = FEPack.meshes.MeshCuboid(1, [0 -1], [0  1], [0 1], numNodesXneg, numNodesYneg, numNodesZneg);
  cell_neg = mesh_neg.domain('volumic');

  % Boundary conditions
  if strcmpi(bc_basis_functions_neg, 'Lagrange')
    init_BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domain('xmin'));
    init_BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domain('xmax'));
  else
    FourierIds = [0 0 0]; FourierIds(1) = floor(numNodesXneg/4);
    init_BCstruct_neg.spB0 = FEPack.spaces.FourierBasis(mesh_neg.domain('xmin'), FourierIds);
    init_BCstruct_neg.spB1 = FEPack.spaces.FourierBasis(mesh_neg.domain('xmax'), FourierIds);
  end

  %% ========= Interior domain
  mesh_int = FEPack.meshes.MeshCuboid(1, BBint, [0 1], [0, 1], numNodesXint, numNodesYint, numNodesZint);
  cell_int  = mesh_int.domain('volumic');
  edgeXminInt = mesh_int.domain('xmin');
  edgeXmaxInt = mesh_int.domain('xmax');
  edge0yInt = mesh_int.domain('ymin');
  edge1yInt = mesh_int.domain('ymax');
  edge0zInt = mesh_int.domain('zmin');
  edge1zInt = mesh_int.domain('zmax');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Finite Element matrices
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  u = FEPack.pdes.PDEObject; 
  v = dual(u);

  %% (+) half-guide
  op_pos.gradu_gradv = (op_pos.funP * (op_pos.Tmat' * grad3(u))) * (op_pos.Tmat' * grad3(v));
  op_pos.gradu_vectv = (op_pos.funP * (op_pos.Tmat' * grad3(u))) * (op_pos.vect * v);
  op_pos.vectu_gradv = (op_pos.funP * (op_pos.vect * u)) * (op_pos.Tmat' * grad3(v));
  op_pos.vectu_vectv = (op_pos.funP * (op_pos.vect * u)) * (op_pos.vect * v);
  op_pos.funQ_u_v    = (op_pos.funQ * u) * v;
  op_pos.funR_u_v    = (op_pos.funR * u) * v;

  fprintf('(+) half-guide FE matrices computation\n');
  mat_pos_gradu_gradv = FEPack.pdes.Form.intg(cell_pos, op_pos.gradu_gradv);
  mat_pos_gradu_vectv = FEPack.pdes.Form.intg(cell_pos, op_pos.gradu_vectv);
  mat_pos_vectu_gradv = FEPack.pdes.Form.intg(cell_pos, op_pos.vectu_gradv);
  mat_pos_vectu_vectv = FEPack.pdes.Form.intg(cell_pos, op_pos.vectu_vectv);
  mat_pos_funQ_u_v    = FEPack.pdes.Form.intg(cell_pos, op_pos.funQ_u_v);
  mat_pos_funR_u_v    = FEPack.pdes.Form.intg(cell_pos, op_pos.funR_u_v);
  fprintf('(+) half-guide FE matrices computation done.\n');

  %% (-) half-guide
  op_neg.gradu_gradv = (op_neg.funP * (op_neg.Tmat' * grad3(u))) * (op_neg.Tmat' * grad3(v));
  op_neg.gradu_vectv = (op_neg.funP * (op_neg.Tmat' * grad3(u))) * (op_neg.vect * v);
  op_neg.vectu_gradv = (op_neg.funP * (op_neg.vect * u)) * (op_neg.Tmat' * grad3(v));
  op_neg.vectu_vectv = (op_neg.funP * (op_neg.vect * u)) * (op_neg.vect * v);
  op_neg.funQ_u_v    = (op_neg.funQ * u) * v;
  op_neg.funR_u_v    = (op_neg.funR * u) * v;

  fprintf('(-) half-guide FE matrices computation\n');
  mat_neg_gradu_gradv = FEPack.pdes.Form.intg(cell_neg, op_neg.gradu_gradv);
  mat_neg_gradu_vectv = FEPack.pdes.Form.intg(cell_neg, op_neg.gradu_vectv);
  mat_neg_vectu_gradv = FEPack.pdes.Form.intg(cell_neg, op_neg.vectu_gradv);
  mat_neg_vectu_vectv = FEPack.pdes.Form.intg(cell_neg, op_neg.vectu_vectv);
  mat_neg_funQ_u_v    = FEPack.pdes.Form.intg(cell_neg, op_neg.funQ_u_v);
  mat_neg_funR_u_v    = FEPack.pdes.Form.intg(cell_neg, op_neg.funR_u_v);
  fprintf('(-) half-guide FE matrices computation done.\n');

  %% Interior domain
  op_int.gradu_gradv = (op_int.funP * (op_int.Tmat' * grad3(u))) * (op_int.Tmat' * grad3(v));
  op_int.gradu_vectv = (op_int.funP * (op_int.Tmat' * grad3(u))) * (op_int.vect * v);
  op_int.vectu_gradv = (op_int.funP * (op_int.vect * u)) * (op_int.Tmat' * grad3(v));
  op_int.vectu_vectv = (op_int.funP * (op_int.vect * u)) * (op_int.vect * v);
  op_int.funQ_u_v    = (op_int.funQ * u) * v;
  op_int.funR_u_v    = (op_int.funR * u) * v;

  fprintf('Interior domain FE matrices computation\n');
  mat_int_gradu_gradv = FEPack.pdes.Form.intg(cell_int, op_int.gradu_gradv);
  mat_int_gradu_vectv = FEPack.pdes.Form.intg(cell_int, op_int.gradu_vectv);
  mat_int_vectu_gradv = FEPack.pdes.Form.intg(cell_int, op_int.vectu_gradv);
  mat_int_vectu_vectv = FEPack.pdes.Form.intg(cell_int, op_int.vectu_vectv);
  mat_int_funQ_u_v    = FEPack.pdes.Form.intg(cell_int, op_int.funQ_u_v);
  mat_int_funR_u_v    = FEPack.pdes.Form.intg(cell_int, op_int.funR_u_v);
  fprintf('Interior domain FE matrices computation done.\n');

  ecs_int = (((u|edge0yInt) - (u|edge1yInt)) == 0.0) & (((u|edge0zInt) - (u|edge1zInt)) == 0.0);
  ecs_int.applyEcs;
  PPint = ecs_int.P;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Edge states
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Variables
  opts.computeSol = false;
  opts.solBasis = false;
  opts.verbose = 0;
  opts.stop_if_spectral_radius_is_one = true;
  err_prefix = '[propagationOperators.m][spectral_radius_is_one]';
  num_err = length(err_prefix);

  % Output initialization
  numD = dispmesh.numPoints;
  val = zeros(numD, 1);

  % Loop
  if (parallel_use)

    if (nargin < 21)
      error('parallel_use has been set to true; you have to specify the number of cores through parallel_numcores. Use the command feature(''numcores'') to see the number of available cores.');
    end

    % Parallel version
    pool = gcp('nocreate'); % Get current parallel pool
    if isempty(pool)
      parpool('local', parallel_numcores); % Start a parallel pool if none exists
    end
    fprintf('Using parfor (parallel execution)\n');

    kpars   = dispmesh.points(:, 1);
    lambdas = dispmesh.points(:, 2);

    parfor idI = 1:numD
      fprintf('%d sur %d\n', idI, numD);

      kpar   = kpars(idI);
      lambda = lambdas(idI); 

      try

        % ========= (+) half-guide
        tic;
        % FE matrix
        AApos = mat_pos_gradu_gradv...
              + 1i * kpar     * mat_pos_vectu_gradv...
              - 1i * kpar     * mat_pos_gradu_vectv...
              + (kpar * kpar) * mat_pos_vectu_vectv...
              + mat_pos_funQ_u_v...
              - lambda * mat_pos_funR_u_v;

        % Boundary condition
        BCstruct_pos = init_BCstruct_pos;
        BCstruct_pos.BCu  = 1i*lambda;
        BCstruct_pos.BCdu = 1;
        %
        % BCstruct_pos.BCu  = 1;
        % BCstruct_pos.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_pos, boundary_op_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 1, AApos, BCstruct_pos, 0, opts);
        tps = toc;
        fprintf('%d (+) guide elapsed time: %f seconds.\n', tps);

        % ========= (-) half-guide
        % FE matrix
        tic;
        AAneg = mat_neg_gradu_gradv...
              + 1i * kpar     * mat_neg_vectu_gradv...
              - 1i * kpar     * mat_neg_gradu_vectv...
              + (kpar * kpar) * mat_neg_vectu_vectv...
              + mat_neg_funQ_u_v...
              - lambda * mat_neg_funR_u_v;

        % Boundary condition
        BCstruct_neg = init_BCstruct_neg;
        BCstruct_neg.BCu  = -1i*lambda;
        BCstruct_neg.BCdu = 1;
        %
        % BCstruct_neg.BCu  = 1;
        % BCstruct_neg.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_neg, boundary_op_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 1, AAneg, BCstruct_neg, 0, opts);
        tps = toc;
        fprintf('%d (-) guide elapsed time: %f seconds.\n', idI, tps);

        % ========= Interior problem
        tic;
        AAint = mat_int_gradu_gradv +...
              +  1i * kpar    * mat_int_vectu_gradv...
              -  1i * kpar    * mat_int_gradu_vectv...
              + (kpar * kpar) * mat_int_vectu_vectv...
              + mat_int_funQ_u_v...
              - lambda * mat_int_funR_u_v;

        % Coefficients associated to boundary condition
        Nbpos = size(boundary_op_pos, 1);
        Nbneg = size(boundary_op_neg, 1);

        BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * boundary_op_pos)...
                    \ (boundary_op_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
        BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * boundary_op_neg)...
                    \ (boundary_op_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

        SSpos = FEPack.pdes.Form.intg_TU_V(edgeXmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
        SSneg = FEPack.pdes.Form.intg_TU_V(edgeXminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

        AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
        tps = toc;
        fprintf('%d interior domain elapsed time: %f seconds.\n', idI, tps);

        % Compute infinity norm of resolvent
        tic;
        % val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));
        % val(idI) = det(AAint0); % SURTOUT PAS !
        val(idI) = condest(AAint0);
        tps = toc;
        fprintf('%d interior domain resolvent elapsed time: %f seconds.\n', idI, tps);

      catch err_obj
        
        % val(idI) = 0;
        if (length(err_obj.message) >= num_err) &&...
          (strcmp(err_obj.message(1:num_err), err_prefix))
          % Some eigenvalues of the Riccati operator are on the unit circle
          val(idI) = 0;
        else
          % Unknown error
          rethrow(err_obj);
        end
      
      end
    end 

  else 

    kpars   = dispmesh.points(:, 1);
    lambdas = dispmesh.points(:, 2);

    for idI = 1:numD
      fprintf('%d sur %d\n', idI, numD);

      kpar   = kpars(idI);
      lambda = lambdas(idI); 

      try

        % ========= (+) half-guide
        tic;
        % FE matrix
        AApos = mat_pos_gradu_gradv...
              + 1i * kpar     * mat_pos_vectu_gradv...
              - 1i * kpar     * mat_pos_gradu_vectv...
              + (kpar * kpar) * mat_pos_vectu_vectv...
              + mat_pos_funQ_u_v...
              - lambda * mat_pos_funR_u_v;

        % Boundary condition
        BCstruct_pos = init_BCstruct_pos;
        BCstruct_pos.BCu  = 1i*lambda;
        BCstruct_pos.BCdu = 1;
        %
        % BCstruct_pos.BCu  = 1;
        % BCstruct_pos.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_pos, boundary_op_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 1, AApos, BCstruct_pos, 0, opts);
        tps = toc;
        fprintf('%d (+) guide elapsed time: %f seconds.\n', tps);

        % ========= (-) half-guide
        % FE matrix
        tic;
        AAneg = mat_neg_gradu_gradv...
              + 1i * kpar     * mat_neg_vectu_gradv...
              - 1i * kpar     * mat_neg_gradu_vectv...
              + (kpar * kpar) * mat_neg_vectu_vectv...
              + mat_neg_funQ_u_v...
              - lambda * mat_neg_funR_u_v;

        % Boundary condition
        BCstruct_neg = init_BCstruct_neg;
        BCstruct_neg.BCu  = -1i*lambda;
        BCstruct_neg.BCdu = 1;
        %
        % BCstruct_neg.BCu  = 1;
        % BCstruct_neg.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_neg, boundary_op_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 1, AAneg, BCstruct_neg, 0, opts);
        tps = toc;
        fprintf('%d (-) guide elapsed time: %f seconds.\n', idI, tps);

        % ========= Interior problem
        tic;
        AAint = mat_int_gradu_gradv +...
              +  1i * kpar    * mat_int_vectu_gradv...
              -  1i * kpar    * mat_int_gradu_vectv...
              + (kpar * kpar) * mat_int_vectu_vectv...
              + mat_int_funQ_u_v...
              - lambda * mat_int_funR_u_v;

        % Coefficients associated to boundary condition
        Nbpos = size(boundary_op_pos, 1);
        Nbneg = size(boundary_op_neg, 1);

        BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * boundary_op_pos)...
                    \ (boundary_op_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
        BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * boundary_op_neg)...
                    \ (boundary_op_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

        SSpos = FEPack.pdes.Form.intg_TU_V(edgeXmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
        SSneg = FEPack.pdes.Form.intg_TU_V(edgeXminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

        AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
        tps = toc;
        fprintf('%d interior domain elapsed time: %f seconds.\n', idI, tps);

        % Compute infinity norm of resolvent
        tic;
        % val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));
        % val(idI) = det(AAint0); % SURTOUT PAS !
        val(idI) = condest(AAint0);
        tps = toc;
        fprintf('%d interior domain resolvent elapsed time: %f seconds.\n', idI, tps);

      catch err_obj
        
        % val(idI) = 0;
        if (length(err_obj.message) >= num_err) &&...
          (strcmp(err_obj.message(1:num_err), err_prefix))
          % Some eigenvalues of the Riccati operator are on the unit circle
          val(idI) = 0;
        else
          % Unknown error
          rethrow(err_obj);
        end
      
      end
    end 

  end

end