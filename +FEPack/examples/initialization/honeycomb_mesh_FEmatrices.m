import FEPack.*

%% Honeycomb lattice
HcObj = tools.HoneycombObject('none', 'none');

% centers = [-1/sqrt(3), 1/sqrt(3); 0, 0];
% ampsV   = [0;  0];
% ampsW   = [10; 10];
% rads    = [0.2; 0.2];
% HcObj.V = @(x) FEPack.tools.FUN_atomic_potential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsV, rads);
% HcObj.W = @(x) 1 + FEPack.tools.FUN_atomic_potential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsW, rads);

HcObj.V = @(x) 0 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
               0 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
               0 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

HcObj.W = @(x) 0 * sin(x(:, 1:2) *  HcObj.dualVec1) +...
               0 * sin(x(:, 1:2) *  HcObj.dualVec2) +...
               0 * sin(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

%% Edge
HcObj.edge.a1 = 1;
HcObj.edge.b1 = 1;
HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);
pbinputs.delta = 1;

%% Meshes and boundary conditions
% # nodes per X/Y edge for periodicity cells
pbinputs.numNodesXpos = 32;
pbinputs.numNodesYpos = pbinputs.numNodesXpos; 
pbinputs.numNodesXneg = pbinputs.numNodesXpos;
pbinputs.numNodesYneg = pbinputs.numNodesXpos;
pbinputs.numNodesXint = pbinputs.numNodesXpos;
pbinputs.numNodesYint = pbinputs.numNodesXpos;

pbinputs.numCellsPos = 2;
pbinputs.numCellsNeg = 2;

% Basis functions for boundary conditions ("Lagrange" or "Fourier")
pbinputs.bc_basis_functions_pos = "Lagrange";
pbinputs.bc_basis_functions_neg = "Lagrange";
pbinputs.bc_basis_functions_int = "Lagrange";

[mshPos, mshNeg, mshInt, op_pos, op_neg, op_int, pbinputs] = construct_mesh_and_FEmatrices_honeycomb(pbinputs, HcObj);

%% Auxiliary functions
function [mshPos, mshNeg, mshInt, op_pos, op_neg, op_int, pbinputsNew] = construct_mesh_and_FEmatrices_honeycomb(pbinputs, HcObj)
  
  %% Coefficients
  % Make sure the edge coefficients are coprime integers
  a1 = HcObj.edge.a1;
  b1 = HcObj.edge.b1;
  [G, x, y] = gcd(a1, b1);  % a1*x + b1*y = G
  a2 = -y; b2 = x;

  if (G ~= 1)
    error('a1 et b1 doivent être premiers entre eux.');
  end

  % Edge vectors
  edge_vec1      =  a1 * HcObj.vecPer1  + b1 * HcObj.vecPer2;
  edge_vec2      =  a2 * HcObj.vecPer1  + b2 * HcObj.vecPer2;
  edge_dual_vec1 =  b2 * HcObj.dualVec1 - a2 * HcObj.dualVec2;
  edge_dual_vec2 = -b1 * HcObj.dualVec1 + a1 * HcObj.dualVec2;

  % Transformation matrices
  Rmat = [edge_vec1, edge_vec2];
  % Tmat = eye(2); warning('Pense à remettre Tmat.');
  Tmat = [edge_dual_vec1'; edge_dual_vec2'] / (2*pi); % Inverse of Rmat
  Tvec = Tmat' * [1; 0];

  % Potentials
  Vpot = @(x) HcObj.V((Rmat * x(:, 1:2)')');
  Wpot = @(x) HcObj.W((Rmat * x(:, 1:2)')');

  % (+) half-guide
  op_pos.funQ = @(x) Vpot(x) + pbinputs.delta * Wpot(x);
  op_pos.funP = 1;
  op_pos.funR = 1;

  % (-) half-guide
  op_neg.funQ = @(x) Vpot(x) - pbinputs.delta * Wpot(x);
  op_neg.funP = 1;
  op_neg.funR = 1;

  % Interior domain
  yInt = ceil(1/(2*pi*pbinputs.delta));
  IDb  = [-yInt, yInt];
  op_int.funQ = @(x) Vpot(x) + pbinputs.delta * HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Wpot(x);
  op_int.funP = 1;
  op_int.funR = 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Translate the functions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  funQposTrans = @(x) op_pos.funQ([x(:, 1) + IDb(2), x(:, 2)]);
  funQnegTrans = @(x) op_neg.funQ([x(:, 1) + IDb(1), x(:, 2)]);

  %% Initialization
  numNodesXpos = pbinputs.numNodesXpos;
  numNodesYpos = pbinputs.numNodesYpos;
  numNodesXneg = pbinputs.numNodesXneg;
  numNodesYneg = pbinputs.numNodesYneg;
  numNodesXint = pbinputs.numNodesXint;
  numNodesYint = pbinputs.numNodesYint;

  bc_basis_functions_pos = pbinputs.bc_basis_functions_pos;
  bc_basis_functions_neg = pbinputs.bc_basis_functions_neg;
  bc_basis_functions_int = pbinputs.bc_basis_functions_int;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Meshes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ========= (+) half-guide
  % Mesh
  mshPos = FEPack.meshes.MeshRectangle(1, [0 1], [0  1], numNodesXpos, numNodesYpos);
  cell_pos = mshPos.domain('volumic');
  
  % Boundary conditions
  if strcmpi(bc_basis_functions_pos, 'Lagrange')
    op_pos.BCstruct.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mshPos.domains{4});
    op_pos.BCstruct.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mshPos.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = floor(numNodesXpos/4);
    op_pos.BCstruct.spB0 = FEPack.spaces.FourierBasis(mshPos.domains{4}, FourierIds);
    op_pos.BCstruct.spB1 = FEPack.spaces.FourierBasis(mshPos.domains{3}, FourierIds);
  end

  %% ========= (-) half-guide
  % Mesh
  mshNeg = FEPack.meshes.MeshRectangle(1, [0 1], [0 -1], numNodesXneg, numNodesYneg);
  cell_neg = mshNeg.domain('volumic');
  
  % Boundary conditions
  if strcmpi(bc_basis_functions_neg, 'Lagrange')
    op_neg.BCstruct.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mshNeg.domains{4});
    op_neg.BCstruct.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mshNeg.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = floor(numNodesXneg/4);
    op_neg.BCstruct.spB0 = FEPack.spaces.FourierBasis(mshNeg.domains{4}, FourierIds);
    op_neg.BCstruct.spB1 = FEPack.spaces.FourierBasis(mshNeg.domains{3}, FourierIds);
  end

  %% ========= Interior domain
  numNodesYint = ceil(numNodesYint * abs(IDb(2) - IDb(1)));
  mshInt = FEPack.meshes.MeshRectangle(1, [0 1], IDb, numNodesXint, numNodesYint);
  cell_int = mshInt.domain('volumic');

  if strcmpi(bc_basis_functions_int, 'Lagrange')
    op_int.spBint_neg = FEPack.spaces.PeriodicLagrangeBasis(mshInt.domains{4});
    op_int.spBint_pos = FEPack.spaces.PeriodicLagrangeBasis(mshInt.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = floor(numNodesXint/4);
    op_int.spBint_pos = FEPack.spaces.FourierBasis(mshInt.domains{4}, FourierIds);
    op_int.spBint_neg = FEPack.spaces.FourierBasis(mshInt.domains{3}, FourierIds);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Finite Element matrices
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  u = FEPack.pdes.PDEObject; 
  v = dual(u);

  %% (+) half-guide
  gradu_gradv_pos = (op_pos.funP * (Tmat' * grad2(u))) * (Tmat' * grad2(v));
  gradu_vectv_pos = (op_pos.funP * (Tmat' * grad2(u))) * (Tvec * v);
  vectu_gradv_pos = (op_pos.funP * (Tvec * u)) * (Tmat' * grad2(v));
  vectu_vectv_pos = (op_pos.funP * (Tvec * u)) * (Tvec * v);
  funQ_u_v_pos    = (funQposTrans * u) * v;
  funR_u_v_pos    = (op_pos.funR * u) * v;

  fprintf('(+) half-guide FE matrices computation\n');
  op_pos.mat_gradu_gradv = FEPack.pdes.Form.intg(cell_pos, gradu_gradv_pos);
  op_pos.mat_gradu_vectv = FEPack.pdes.Form.intg(cell_pos, gradu_vectv_pos);
  op_pos.mat_vectu_gradv = FEPack.pdes.Form.intg(cell_pos, vectu_gradv_pos);
  op_pos.mat_vectu_vectv = FEPack.pdes.Form.intg(cell_pos, vectu_vectv_pos);
  op_pos.mat_funQ_u_v    = FEPack.pdes.Form.intg(cell_pos, funQ_u_v_pos);
  op_pos.mat_funR_u_v    = FEPack.pdes.Form.intg(cell_pos, funR_u_v_pos);
  fprintf('BCstructNeg(+) half-guide FE matrices computation done.\n');

  %% (-) half-guide
  gradu_gradv_neg = (op_neg.funP * (Tmat' * grad2(u))) * (Tmat' * grad2(v));
  gradu_vectv_neg = (op_neg.funP * (Tmat' * grad2(u))) * (Tvec * v);
  vectu_gradv_neg = (op_neg.funP * (Tvec * u)) * (Tmat' * grad2(v));
  vectu_vectv_neg = (op_neg.funP * (Tvec * u)) * (Tvec * v);
  funQ_u_v_neg    = (funQnegTrans * u) * v;
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
  gradu_gradv_int = (op_int.funP * (Tmat' * grad2(u))) * (Tmat' * grad2(v));
  gradu_vectv_int = (op_int.funP * (Tmat' * grad2(u))) * (Tvec * v);
  vectu_gradv_int = (op_int.funP * (Tvec * u)) * (Tmat' * grad2(v));
  vectu_vectv_int = (op_int.funP * (Tvec * u)) * (Tvec * v);
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

  %% New pbinputs
  pbinputsNew = pbinputs;
  pbinputsNew.Rmat = Rmat;
  pbinputsNew.Tmat = Tmat;
  pbinputsNew.Tvec = Tvec;
  pbinputsNew.IDb = IDb;

end