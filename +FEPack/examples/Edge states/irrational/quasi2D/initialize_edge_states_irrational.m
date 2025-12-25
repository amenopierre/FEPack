function [...
    meshXYpos, meshXYneg, meshYZ, meshLineZ,...
    init_BCstruct_pos, init_BCstruct_neg...
  ]...
   =...
  initialize_edge_states_irrational(...
    numNodesX, numNodesY, numNodesZ, is_mesh_structured,...
    bc_basis_functions_X, bc_basis_functions_S,...
    op_pos, op_neg, cutvec,...
    folder_name,...
    parallel_use, parallel_numcores...
  )

  %% ========= Meshes
  meshXYpos = FEPack.meshes.MeshRectangle(is_mesh_structured, [0  1], [0  1], numNodesX, numNodesY);
  meshXYneg = FEPack.meshes.MeshRectangle(is_mesh_structured, [0 -1], [0  1], numNodesX, numNodesY);
  meshYZ    = FEPack.meshes.MeshRectangle(is_mesh_structured, [0  1], [0  1], numNodesY, numNodesZ);
  meshLineZ = FEPack.meshes.MeshSegment('uniform', 0, 1, numNodesZ);

  %% ========= Boundary-related basis functions
  % (+) guide
  if strcmpi(bc_basis_functions_X, 'Lagrange')
    init_BCstruct_pos.spBX = FEPack.spaces.PeriodicLagrangeBasis(meshYZ.domain('volumic')); % X = cst
  else
    nBY = max(1, floor(numNodesY/4));
    nBZ = max(1, floor(numNodesZ/4));
    init_BCstruct_pos.spBX = FEPack.spaces.FourierBasis(meshYZ.domain('volumic'), [nBY, nBZ]);
  end

  if strcmpi(bc_basis_functions_S, 'Lagrange')
    init_BCstruct_pos.spBS = FEPack.spaces.PeriodicLagrangeBasis(meshLineZ.domain('volumic')); % S/Z line
  else
    nBS = max(1, floor(numNodesZ/4));
    init_BCstruct_pos.spBS = FEPack.spaces.FourierBasis(meshLineZ.domain('volumic'), nBS);
  end

  % (-) guide
  init_BCstruct_neg = init_BCstruct_pos;

  %% ========= Finite Element matrices

  % (+) guide
  compute_FE_matrices(...
    meshXYpos, meshLineZ,...
    op_pos.funP3D, op_pos.funQ3D, op_pos.funR3D, op_pos.Tmat, op_pos.vect, cutvec,...
    folder_name, 'pos',...
    parallel_use, parallel_numcores);

  % (-) guide
  compute_FE_matrices(...
    meshXYneg, meshLineZ,...
    op_neg.funP3D, op_neg.funQ3D, op_neg.funR3D, op_neg.Tmat, op_neg.vect, cutvec,...
    folder_name, 'neg',...
    parallel_use, parallel_numcores);
  
end

% Auxiliary functions
function compute_FE_matrices(...
  meshXY, meshLineZ,...
  funP3D, funQ3D, funR3D, Tmat, vect, cutvec,...
  folder_name, out_suffix,...
  parallel_use, parallel_numcores)

  u = FEPack.pdes.PDEObject; 
  v = dual(u);
  
  % Integrands for FE matrices
  gradu_gradv = @(funP) (funP * (Tmat' * grad2(u))) * (Tmat' * grad2(v));
  gradu_vectv = @(funP) (funP * (Tmat' * grad2(u))) * (vect * v);
  vectu_gradv = @(funP) (funP * (vect * u)) * (Tmat' * grad2(v));
  vectu_vectv = @(funP) (funP * (vect * u)) * (vect * v);
  funQ_u_v    = @(funQ) (funQ * u) * v;
  funR_u_v    = @(funR) (funR * u) * v;

  % FE matrices
  N_s = meshLineZ.numPoints;
  sVec = meshLineZ.points;
  % N_XY = meshXY.numPoints;
  cellXY = meshXY.domain('volumic');
  
  % fid = fopen([nomdossier, 'elapsed_time', suffix], 'w+');
  fid = fopen('outputs/elapsed_time', 'a+');
  fprintf('Finite elements matrices calculations.\n');

  tic;
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

    parfor idS = 1:N_s-1
      sVar = sVec(idS, 1);
      funPs = coefficient_2D_trace(funP3D, cutvec, sVar);
      funQs = coefficient_2D_trace(funQ3D, cutvec, sVar);
      funRs = coefficient_2D_trace(funR3D, cutvec, sVar);
      
      FEmat = [];
      FEmat.mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, gradu_gradv(funPs));
      FEmat.mat_gradu_vectv = FEPack.pdes.Form.intg(cellXY, gradu_vectv(funPs));
      FEmat.mat_vectu_gradv = FEPack.pdes.Form.intg(cellXY, vectu_gradv(funPs));
      FEmat.mat_vectu_vectv = FEPack.pdes.Form.intg(cellXY, vectu_vectv(funPs));
      FEmat.mat_funQ_u_v    = FEPack.pdes.Form.intg(cellXY, funQ_u_v(funQs));
      FEmat.mat_funR_u_v    = FEPack.pdes.Form.intg(cellXY, funR_u_v(funRs));

      FEPack.tools.parsave([folder_name, 'FEmat_', out_suffix, '_', num2str(idS)], FEmat, true);
    end

  else

    for idS = 1:N_s-1
      fprintf('%d sur %d\n', idS, N_s);

      sVar = sVec(idS, 1);
      funPs = coefficient_2D_trace(funP3D, cutvec, sVar);
      funQs = coefficient_2D_trace(funQ3D, cutvec, sVar);
      funRs = coefficient_2D_trace(funR3D, cutvec, sVar);

      FEmat = [];
      FEmat.mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, gradu_gradv(funPs));
      FEmat.mat_gradu_vectv = FEPack.pdes.Form.intg(cellXY, gradu_vectv(funPs));
      FEmat.mat_vectu_gradv = FEPack.pdes.Form.intg(cellXY, vectu_gradv(funPs));
      FEmat.mat_vectu_vectv = FEPack.pdes.Form.intg(cellXY, vectu_vectv(funPs));
      FEmat.mat_funQ_u_v    = FEPack.pdes.Form.intg(cellXY, funQ_u_v(funQs));
      FEmat.mat_funR_u_v    = FEPack.pdes.Form.intg(cellXY, funR_u_v(funRs));

      FEPack.tools.parsave([folder_name, 'FEmat_', out_suffix, '_', num2str(idS)], FEmat, true);
    end

  end

  tps = toc;
  fprintf(fid, '%0.5e\t', tps);
  fclose(fid);
end

function fun_s = coefficient_2D_trace(fun3D, cutvec, s)

  if isa(fun3D, 'function_handle')
    fun_s = @(x) fun3D([x(:, 1), cutvec(1)*x(:, 2), cutvec(2)*x(:, 2) + s]);
  else
    fun_s = fun3D;
  end

end