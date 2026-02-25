function [mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = initialize_edge_states_rational(...
    op_pos, numNodesXpos, numNodesYpos, bc_basis_functions_pos,...
    op_neg, numNodesXneg, numNodesYneg, bc_basis_functions_neg,...
    op_int, BBint, numNodesXint, numNodesYint, ...
    problem_type...
  )

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Meshes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% ========= (+) half-guide
  % Mesh
  mesh_pos = FEPack.meshes.MeshRectangle(1, [0 1], [0  1], numNodesXpos, numNodesYpos);
  cell_pos = mesh_pos.domain('volumic');
  
  % Boundary conditions
  if strcmpi(bc_basis_functions_pos, 'Lagrange')
    BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{4});
    BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = floor(numNodesXpos/4);
    BCstruct_pos.spB0 = FEPack.spaces.FourierBasis(mesh_pos.domains{4}, FourierIds);
    BCstruct_pos.spB1 = FEPack.spaces.FourierBasis(mesh_pos.domains{3}, FourierIds);
  end

  %% ========= (-) half-guide
  % Mesh
  mesh_neg = FEPack.meshes.MeshRectangle(1, [0 1], [0 -1], numNodesXneg, numNodesYneg);
  cell_neg = mesh_neg.domain('volumic');
  
  % Boundary conditions
  if strcmpi(bc_basis_functions_neg, 'Lagrange')
    BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{4});
    BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = floor(numNodesXneg/4);
    BCstruct_neg.spB0 = FEPack.spaces.FourierBasis(mesh_neg.domains{4}, FourierIds);
    BCstruct_neg.spB1 = FEPack.spaces.FourierBasis(mesh_neg.domains{3}, FourierIds);
  end

  %% ========= Interior domain
  if strcmpi(problem_type, 'interior')
    mesh_int = FEPack.meshes.MeshRectangle(1, [0 1], BBint, numNodesXint, numNodesYint);
    cell_int = mesh_int.domain('volumic');
  else
    mesh_int = [];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Finite Element matrices
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  u = FEPack.pdes.PDEObject; 
  v = dual(u);

  %% (+) half-guide
  gradu_gradv_pos = (op_pos.funP * (op_pos.Tmat' * grad2(u))) * (op_pos.Tmat' * grad2(v));
  gradu_vectv_pos = (op_pos.funP * (op_pos.Tmat' * grad2(u))) * (op_pos.vect * v);
  vectu_gradv_pos = (op_pos.funP * (op_pos.vect * u)) * (op_pos.Tmat' * grad2(v));
  vectu_vectv_pos = (op_pos.funP * (op_pos.vect * u)) * (op_pos.vect * v);
  funQ_u_v_pos    = (op_pos.funQ * u) * v;
  funR_u_v_pos    = (op_pos.funR * u) * v;

  fprintf('(+) half-guide FE matrices computation\n');
  op_pos.mat_gradu_gradv = FEPack.pdes.Form.intg(cell_pos, gradu_gradv_pos);
  op_pos.mat_gradu_vectv = FEPack.pdes.Form.intg(cell_pos, gradu_vectv_pos);
  op_pos.mat_vectu_gradv = FEPack.pdes.Form.intg(cell_pos, vectu_gradv_pos);
  op_pos.mat_vectu_vectv = FEPack.pdes.Form.intg(cell_pos, vectu_vectv_pos);
  op_pos.mat_funQ_u_v    = FEPack.pdes.Form.intg(cell_pos, funQ_u_v_pos);
  op_pos.mat_funR_u_v    = FEPack.pdes.Form.intg(cell_pos, funR_u_v_pos);
  fprintf('(+) half-guide FE matrices computation done.\n');

  %% (-) half-guide
  gradu_gradv_neg = (op_neg.funP * (op_neg.Tmat' * grad2(u))) * (op_neg.Tmat' * grad2(v));
  gradu_vectv_neg = (op_neg.funP * (op_neg.Tmat' * grad2(u))) * (op_neg.vect * v);
  vectu_gradv_neg = (op_neg.funP * (op_neg.vect * u)) * (op_neg.Tmat' * grad2(v));
  vectu_vectv_neg = (op_neg.funP * (op_neg.vect * u)) * (op_neg.vect * v);
  funQ_u_v_neg    = (op_neg.funQ * u) * v;
  funR_u_v_neg    = (op_neg.funR * u) * v;

  fprintf('(-) half-guide FE matrices computation\n');
  op_neg.mat_gradu_gradv = FEPack.pdes.Form.intg(cell_neg, gradu_gradv_neg);
  op_neg.mat_gradu_vectv = FEPack.pdes.Form.intg(cell_neg, gradu_vectv_neg);
  op_neg.mat_vectu_gradv = FEPack.pdes.Form.intg(cell_neg, vectu_gradv_neg);
  op_neg.mat_vectu_vectv = FEPack.pdes.Form.intg(cell_neg, vectu_vectv_neg);
  op_neg.mat_funQ_u_v    = FEPack.pdes.Form.intg(cell_neg, funQ_u_v_neg);
  op_neg.mat_funR_u_v    = FEPack.pdes.Form.intg(cell_neg, funR_u_v_neg);
  fprintf('(-) half-guide FE matrices computation done.\n');

  %% Interior domain
  if strcmpi(problem_type, 'interior')
    gradu_gradv_int = (op_int.funP * (op_int.Tmat' * grad2(u))) * (op_int.Tmat' * grad2(v));
    gradu_vectv_int = (op_int.funP * (op_int.Tmat' * grad2(u))) * (op_int.vect * v);
    vectu_gradv_int = (op_int.funP * (op_int.vect * u)) * (op_int.Tmat' * grad2(v));
    vectu_vectv_int = (op_int.funP * (op_int.vect * u)) * (op_int.vect * v);
    funQ_u_v_int    = (op_int.funQ * u) * v;
    funR_u_v_int    = (op_int.funR * u) * v;

    fprintf('Interior domain FE matrices computation\n');
    op_int.mat_gradu_gradv = FEPack.pdes.Form.intg(cell_int, gradu_gradv_int);
    op_int.mat_gradu_vectv = FEPack.pdes.Form.intg(cell_int, gradu_vectv_int);
    op_int.mat_vectu_gradv = FEPack.pdes.Form.intg(cell_int, vectu_gradv_int);
    op_int.mat_vectu_vectv = FEPack.pdes.Form.intg(cell_int, vectu_vectv_int);
    op_int.mat_funQ_u_v    = FEPack.pdes.Form.intg(cell_int, funQ_u_v_int);
    op_int.mat_funR_u_v    = FEPack.pdes.Form.intg(cell_int, funR_u_v_int);
    fprintf('Interior domain FE matrices computation done.\n');
  else
    op_int.mat_gradu_gradv = [];
    op_int.mat_gradu_vectv = [];
    op_int.mat_vectu_gradv = [];
    op_int.mat_vectu_vectv = [];
    op_int.mat_funQ_u_v    = [];
    op_int.mat_funR_u_v    = [];
  end
end