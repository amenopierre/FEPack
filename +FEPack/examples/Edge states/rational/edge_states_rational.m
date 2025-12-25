function [val, eigvals] = edge_states_rational(...
  dispersion_points,...
  op_pos, mesh_pos, init_BCstruct_pos,...
  op_neg, mesh_neg, init_BCstruct_neg,...
  op_int, mesh_int, problem_type,...
  parallel_use, parallel_numcores, plot_coefficients, stop_at_nans)
  
  %% [mesh, val] = EDGE_STATES_RATIONAL()
  % 

  if (nargin < 14)
    stop_at_nans = false;
  end

  %% Plot coefficients
  if (nargin >= 13 && plot_coefficients)
    % Plot coefficient
    numX = 2;
    numY = 2;

    % if isa(op_int.funP, 'function_handle')
    %   funP = op_int.funP;
    % else
    %   funP = @(x) ones(size(x, 1), 1) * op_int.funP;
    % end

    if isa(op_int.funQ, 'function_handle')
      funQ = op_int.funQ;
    else
      funQ = @(x) ones(size(x, 1), 1) * op_int.funQ;
    end
    
    figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    for idX = -numX:numX-1
      for idY = -numY:numY-1
        X = mesh_pos.points(:, 1) + idX;
        Y = mesh_pos.points(:, 2) + idY;

        % % P
        % subplot(1, 2, 1);
        % trisurf(mesh_pos.triangles, X, Y, funP([X, Y]));
        % hold on;
        % shading interp;
        % view(2);
        % set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
        % colorbar('TickLabelInterpreter', 'latex');
        % xlabel('$P$');

        % Q
        % subplot(1, 2, 2);
        trisurf(mesh_pos.triangles, X, Y, funQ([X, Y]));
        hold on;
        shading interp;
        view(2);
        set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
        colorbar('TickLabelInterpreter', 'latex');
        xlabel('$Q$');
      end
    end

    pause;
  end

  %% Finite Element matrices
  mat_pos_gradu_gradv = op_pos.mat_gradu_gradv;
  mat_pos_gradu_vectv = op_pos.mat_gradu_vectv;
  mat_pos_vectu_gradv = op_pos.mat_vectu_gradv;
  mat_pos_vectu_vectv = op_pos.mat_vectu_vectv;
  mat_pos_funQ_u_v    = op_pos.mat_funQ_u_v;
  mat_pos_funR_u_v    = op_pos.mat_funR_u_v;

  mat_neg_gradu_gradv = op_neg.mat_gradu_gradv;
  mat_neg_gradu_vectv = op_neg.mat_gradu_vectv;
  mat_neg_vectu_gradv = op_neg.mat_vectu_gradv;
  mat_neg_vectu_vectv = op_neg.mat_vectu_vectv;
  mat_neg_funQ_u_v    = op_neg.mat_funQ_u_v;
  mat_neg_funR_u_v    = op_neg.mat_funR_u_v;

  if strcmpi(problem_type, 'interior')
    mat_int_gradu_gradv = op_int.mat_gradu_gradv;
    mat_int_gradu_vectv = op_int.mat_gradu_vectv;
    mat_int_vectu_gradv = op_int.mat_vectu_gradv;
    mat_int_vectu_vectv = op_int.mat_vectu_vectv;
    mat_int_funQ_u_v    = op_int.mat_funQ_u_v;
    mat_int_funR_u_v    = op_int.mat_funR_u_v;

    % Essential conditions for interior problem
    edge0xInt   = mesh_int.domain('xmin');
    edge1xInt   = mesh_int.domain('xmax');
    edgeYminInt = mesh_int.domain('ymin');
    edgeYmaxInt = mesh_int.domain('ymax');
    u = FEPack.pdes.PDEObject; 
    ecs_int = (((u|edge0xInt) - (u|edge1xInt)) == 0.0);
    ecs_int.applyEcs;
    PPint = ecs_int.P;
  else
    % This is apparently important when using parallel computing (parfor)
    mat_int_gradu_gradv = [];
    mat_int_gradu_vectv = [];
    mat_int_vectu_gradv = [];
    mat_int_vectu_vectv = [];
    mat_int_funQ_u_v    = [];
    mat_int_funR_u_v    = [];
    edgeYmaxInt = [];
    edgeYminInt = [];
    PPint = [];
  end

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
  numD = size(dispersion_points, 1);
  val = zeros(numD, 1);
  eigvals = cell(numD, 1);

  % Loop
  if (parallel_use)

    if (nargin < 12)
      error('parallel_use has been set to true; you have to specify the number of cores through parallel_numcores. Use the command feature(''numcores'') to see the number of available cores.');
    end

    % Parallel version
    pool = gcp('nocreate'); % Get current parallel pool
    if isempty(pool)
      parpool('local', parallel_numcores); % Start a parallel pool if none exists
    end
    fprintf('Using parfor (parallel execution)\n');

    kpars   = dispersion_points(:, 1);
    lambdas = dispersion_points(:, 2);

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

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_pos, boundary_op_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 2, AApos, BCstruct_pos, 0, opts);
        tps = toc;
        fprintf('(+) guide elapsed time: %f seconds.\n', tps);

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

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_neg, boundary_op_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 2, AAneg, BCstruct_neg, 0, opts);
        tps = toc;
        fprintf('(-) guide elapsed time: %f seconds.\n', idI, tps);
        
        if ~(stop_at_nans)
          if strcmpi(problem_type, 'interior')

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

            SSpos = FEPack.pdes.Form.intg_TU_V(edgeYmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
            SSneg = FEPack.pdes.Form.intg_TU_V(edgeYminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

            AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
            tps = toc;
            fprintf('interior domain elapsed time: %f seconds.\n', idI, tps);

            % Compute infinity norm of resolvent
            tic;
            val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));
            % val(idI) = det(AAint0); % SURTOUT PAS !
            % val(idI) = condest(AAint0);
            eigvals{idI} = eig(AAint0).';
            tps = toc;
            fprintf('interior domain resolvent elapsed time: %f seconds.\n', idI, tps);

          elseif strcmpi(problem_type, 'interface')

            % ========= Interface problem
            tic;
            val(idI) = cond(boundary_op_pos + boundary_op_neg);
            eigvals{idI} = eig(boundary_op_pos + boundary_op_neg);
            tps = toc;
            fprintf('interface problem elapsed time: %f seconds.\n', idI, tps);
            
          else
            error(['Unrecognized value ''', problem_type, ''' for problem_type. Possible values are ''interior'' and ''interface''.']);
          end
        else
          eigvals{idI} = NaN;
        end

      catch err_obj
        
        % val(idI) = 0;
        if (length(err_obj.message) >= num_err) &&...
          (strcmp(err_obj.message(1:num_err), err_prefix))
          % Some eigenvalues of the Riccati operator are on the unit circle
          val(idI) = NaN;
          eigvals{idI} = NaN;
        else
          % Unknown error
          rethrow(err_obj);
        end
      
      end
    end 

  else 

    kpars   = dispersion_points(:, 1);
    lambdas = dispersion_points(:, 2);

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
        BCstruct_pos.BCu  = 1; % 1i*lambda;
        BCstruct_pos.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_pos, boundary_op_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 2, AApos, BCstruct_pos, 0, opts);
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
        BCstruct_neg.BCu  = 1;%-1i*lambda;
        BCstruct_neg.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_neg, boundary_op_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 2, AAneg, BCstruct_neg, 0, opts);
        tps = toc;
        fprintf('%d (-) guide elapsed time: %f seconds.\n', idI, tps);

        if ~(stop_at_nans)
          if strcmpi(problem_type, 'interior')
            
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

            SSpos = FEPack.pdes.Form.intg_TU_V(edgeYmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
            SSneg = FEPack.pdes.Form.intg_TU_V(edgeYminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

            AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
            tps = toc;
            fprintf('%d interior domain elapsed time: %f seconds.\n', idI, tps);

            % Compute infinity norm of resolvent
            tic;
            % val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));
            % val(idI) = det(AAint0); % SURTOUT PAS !
            val(idI) = condest(AAint0);
            eigvals{idI} = eigs(AAint0, size(AAint0, 1)).';
            tps = toc;
            fprintf('%d interior domain resolvent elapsed time: %f seconds.\n', idI, tps);

          elseif strcmpi(problem_type, 'interface')

            % ========= Interface problem
            tic;
            val(idI) = cond(boundary_op_pos + boundary_op_neg);
            eigvals{idI} = eig(boundary_op_pos + boundary_op_neg);
            tps = toc;
            fprintf('%d interface problem elapsed time: %f seconds.\n', idI, tps);
          
          else
            error(['Unrecognized value ''', problem_type, ''' for problem_type. Possible values are ''interior'' and ''interface''.']);
          end
        else
          eigvals{idI} = NaN;
        end

      catch err_obj
        
        % val(idI) = 0;
        if (length(err_obj.message) >= num_err) &&...
          (strcmp(err_obj.message(1:num_err), err_prefix))
          % Some eigenvalues of the Riccati operator are on the unit circle
          val(idI) = NaN;
          eigvals{idI} = NaN;
        else
          % Unknown error
          rethrow(err_obj);
        end
      
      end
    end 

  end

end