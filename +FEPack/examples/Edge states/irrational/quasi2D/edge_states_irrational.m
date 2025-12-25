function [dispmesh, val] = edge_states_irrational(...
  BBkpar, numNodesK,...
  BBlambda, numNodesLambda,...
  meshXYpos, meshXYneg, meshYZ, meshLineZ, quadratureRatio,...
  cutvec,...
  init_BCstruct_pos, init_BCstruct_neg,...
  folder_name,...
  parallel_use, parallel_numcores)

  %% Mesh associated to parallel quasi-momentum / spectral parameter
  dispmesh = FEPack.meshes.MeshRectangle(0, BBkpar, BBlambda, numNodesK, numNodesLambda);

  %%

  %% Compute edge states
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

    if (nargin < 18)
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
        FEcos = struct('gradu_gradv',  1,...
                       'gradu_vectv',  1i * kpar,...
                       'vectu_gradv', -1i * kpar,...
                       'vectu_vectv',  kpar * kpar,...
                       'funQ_u_v',     1,...
                       'funR_u_v',    -lambda);
        
        % Boundary condition
        BCstruct_pos = init_BCstruct_pos;
        BCstruct_pos.BCu  = 1i;%*lambda;
        BCstruct_pos.BCdu = 1;

        % Half-guide problem
        [~, ~, ~, boundary_op_pos] = quasi2DHalfGuide(meshXYpos, meshYZ, meshLineZ, +1, FEcos, cutvec, BCstruct_pos, quadratureRatio, 'pos', folder_name, opts);
        tps = toc;
        fprintf('%d (+) guide elapsed time: %f seconds.\n', idI, tps);

        % ========= (-) half-guide
        tic;
        FEcos = struct('gradu_gradv',  1,...
                       'gradu_vectv',  1i * kpar,...
                       'vectu_gradv', -1i * kpar,...
                       'vectu_vectv',  kpar * kpar,...
                       'funQ_u_v',     1,...
                       'funR_u_v',    -lambda);
        
        % Boundary condition
        BCstruct_neg = init_BCstruct_neg;
        BCstruct_neg.BCu  = -1i*lambda;
        BCstruct_neg.BCdu = 1;

        % Half-guide problem
        [~, ~, ~, boundary_op_neg] = quasi2DHalfGuide(meshXYneg, meshYZ, meshLineZ, +1, FEcos, cutvec, BCstruct_neg, quadratureRatio, 'neg', folder_name, opts);
        tps = toc;
        fprintf('%d (-) guide elapsed time: %f seconds.\n', idI, tps);

        % ========= Interior problem
        tic;
        val(idI) = cond(boundary_op_pos + boundary_op_neg);
        % AAint = mat_int_gradu_gradv +...
        %       +  1i * kpar    * mat_int_vectu_gradv...
        %       -  1i * kpar    * mat_int_gradu_vectv...
        %       + (kpar * kpar) * mat_int_vectu_vectv...
        %       + mat_int_funQ_u_v...
        %       - lambda * mat_int_funR_u_v;

        % % Coefficients associated to boundary condition
        % Nbpos = size(boundary_op_pos, 1);
        % Nbneg = size(boundary_op_neg, 1);

        % BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * boundary_op_pos)...
        %             \ (boundary_op_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
        % BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * boundary_op_neg)...
        %             \ (boundary_op_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

        % SSpos = FEPack.pdes.Form.intg_TU_V(edgeYmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
        % SSneg = FEPack.pdes.Form.intg_TU_V(edgeYminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

        % AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
        % tps = toc;
        % fprintf('%d interior domain elapsed time: %f seconds.\n', idI, tps);

        % % Compute infinity norm of resolvent
        % tic;
        % % val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));
        % % val(idI) = det(AAint0); % SURTOUT PAS !
        % val(idI) = condest(AAint0);
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
        FEcos = struct('gradu_gradv',  1,...
                       'gradu_vectv',  1i * kpar,...
                       'vectu_gradv', -1i * kpar,...
                       'vectu_vectv',  kpar * kpar,...
                       'funQ_u_v',     1,...
                       'funR_u_v',    -lambda);
        
        % Boundary condition
        BCstruct_pos = init_BCstruct_pos;
        BCstruct_pos.BCu  = 1i*lambda;
        BCstruct_pos.BCdu = 1;

        % Half-guide problem
        [~, ~, ~, boundary_op_pos] = quasi2DHalfGuide(meshXYpos, meshYZ, meshLineZ, +1, FEcos, cutvec, BCstruct_pos, quadratureRatio, 'pos', folder_name, opts);
        tps = toc;
        fprintf('%d (+) guide elapsed time: %f seconds.\n', idI, tps);

        % ========= (-) half-guide
        tic;
        FEcos = struct('gradu_gradv',  1,...
                       'gradu_vectv',  1i * kpar,...
                       'vectu_gradv', -1i * kpar,...
                       'vectu_vectv',  kpar * kpar,...
                       'funQ_u_v',     1,...
                       'funR_u_v',    -lambda);
        
        % Boundary condition
        BCstruct_neg = init_BCstruct_neg;
        BCstruct_neg.BCu  = -1i*lambda;
        BCstruct_neg.BCdu = 1;

        % Half-guide problem
        [~, ~, ~, boundary_op_neg] = quasi2DHalfGuide(meshXYneg, meshYZ, meshLineZ, +1, FEcos, cutvec, BCstruct_neg, quadratureRatio, 'neg', folder_name, opts);
        tps = toc;
        fprintf('%d (-) guide elapsed time: %f seconds.\n', idI, tps);

        % ========= Interior problem
        tic;
        val(idI) = cond(boundary_op_pos + boundary_op_neg);
        % AAint = mat_int_gradu_gradv +...
        %       +  1i * kpar    * mat_int_vectu_gradv...
        %       -  1i * kpar    * mat_int_gradu_vectv...
        %       + (kpar * kpar) * mat_int_vectu_vectv...
        %       + mat_int_funQ_u_v...
        %       - lambda * mat_int_funR_u_v;

        % % Coefficients associated to boundary condition
        % Nbpos = size(boundary_op_pos, 1);
        % Nbneg = size(boundary_op_neg, 1);

        % BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * boundary_op_pos)...
        %             \ (boundary_op_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
        % BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * boundary_op_neg)...
        %             \ (boundary_op_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

        % SSpos = FEPack.pdes.Form.intg_TU_V(edgeYmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
        % SSneg = FEPack.pdes.Form.intg_TU_V(edgeYminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

        % AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
        % tps = toc;
        % fprintf('%d interior domain elapsed time: %f seconds.\n', idI, tps);

        % % Compute infinity norm of resolvent
        % tic;
        % % val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));
        % % val(idI) = det(AAint0); % SURTOUT PAS !
        % val(idI) = condest(AAint0);
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