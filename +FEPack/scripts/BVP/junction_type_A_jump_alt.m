clear; clc;

%% NOTE TO MYSELF: LOOP OVER POS AND NEG!!!
% Preliminary
import FEPack.*
% profile ON

%% Problem-related variables
opts.omega = 8 + 0.5i;
opts.verbose = 0;
opts.computeSol = true;
opts.solBasis = true;

% Periods along the interface
S.pos.periodY = 1;
S.neg.periodY = 0.5*sqrt(2);

% Coefficients
EX = [1; 0];
EY.pos = [0; S.pos.periodY];
EY.neg = [0; S.neg.periodY];
C  = [0.5, 0.5];
S.pos.mu  = @(x) 0.5 + tools.FUN_per_cutoff_circle(x, EX, EY.pos, C, [-0.3, 0.3]);
S.pos.rho = @(x) 0.5 + tools.FUN_per_cutoff_cuboid(x, EX, EY.neg, C, [-0.3, 0.3], [-0.3, 0.3]);

S.neg.mu  = @(x) 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2) / S.neg.periodY);
S.neg.rho = @(x) 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2) / S.neg.periodY);

% Source term
source = @(x) tools.cutoff(sqrt(x(:, 1).^2 + x(:, 2).^2), -0.75, 0.75);

% numbers of cells
S.pos.numCellsX = 2; 
S.pos.numCellsY = 2;
S.pos.numFBpts =  50;

S.neg.numCellsX = 2; 
S.neg.numCellsY = 2;
S.neg.numFBpts =  50; 

S.pos.numFBpts = floor(S.pos.numFBpts * 2*pi / S.pos.periodY);
S.neg.numFBpts = floor(S.neg.numFBpts * 2*pi / S.neg.periodY);

%% Positive vs Negative
sides.name = {'pos', 'neg'};
sides.sign = {+1, -1};

%% Mesh
struct_mesh = 1;
N = 32;
S.pos.numNodesX = N; 
S.pos.numNodesY = floor(S.pos.periodY * N);
S.neg.numNodesX = N; 
S.neg.numNodesY = floor(S.neg.periodY * N); 

for Iside = 1:2
  sn = sides.name{Iside};
  ss = sides.sign{Iside};
  S.(sn).mesh = meshes.MeshRectangle(struct_mesh, ss*EX, EY.(sn), S.(sn).numNodesX, S.(sn).numNodesY);
end

% %% Basis functions associated to the interface
% type_basis = 'Lagrange'; % 'sine', 'Lagrange'

% if strcmpi(type_basis, 'sine')

%   int_basis_functions = @(x, n) sqrt(1 / Lint) * sin(pi * (x(:, 1) + Lint) * n / (2 * Lint)) .* (abs(x(:, 1)) <= Lint);
%   numBasisInt = 20;

% elseif strcmpi(type_basis, 'Lagrange')

%   % error('Il n''y a pas de bug. Je veux juste que tu vérifies le nombre de noeuds (pas trop grand).');
%   int_basis_functions = @(x, n) ((x(:, 1) - mesh2Dint.points(n, 1).') ./ (mesh2Dint.points(n+1, 1).' - mesh2Dint.points(n, 1).'))...
%                                   .* (x(:, 1) >= mesh2Dint.points(n, 1).' & x(:, 1) < mesh2Dint.points(n+1, 1).')...
%                                + ((mesh2Dint.points(n+2, 1).' - x(:, 1)) ./ (mesh2Dint.points(n+2, 1).' - mesh2Dint.points(n+1, 1).'))...
%                                   .* (x(:, 1) >= mesh2Dint.points(n+1, 1).' & x(:, 1) < mesh2Dint.points(n+2, 1).');
%   numBasisInt = mesh2Dint.numPoints - 2;

% else

%   error(['Type de fonction de base', type_basis, 'non reconnu.']);

% end

% spBint = spaces.SpectralBasis(mesh2Dint.domain('volumic'), int_basis_functions, numBasisInt);
% spBint.computeBasisMatrices;


%% Plot coefficients?
plot_coefficients = false;

if (plot_coefficients) 
  set(groot,'defaultAxesTickLabelInterpreter','latex'); %#ok
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  mu  = @(x) (x(:, 1) >= 0) .*  S.pos.mu(x) + (x(:, 1) <  0) .* S.neg.mu(x);
  rho = @(x) (x(:, 1) >= 0) .* S.pos.rho(x) + (x(:, 1) <  0) .* S.neg.rho(x);

  for Iside = 1:2
    sn = sides.name{Iside};
    ss = sides.sign{Iside};

    for idX = 1:S.(sn).numCellsX
      for idY = 1:(2*S.(sn).numCellsY)
        X = S.(sn).mesh.points(:, 1) + ss*(idX - 1);
        Y = S.(sn).mesh.points(:, 2) + (idY - S.(sn).numCellsY - 1)*S.(sn).periodY;

        figure (1);
        trisurf(S.(sn).mesh.triangles, X, Y, mu([X, Y]));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

        figure (2);
        trisurf(S.(sn).mesh.triangles, X, Y, rho([X, Y]));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

        figure (3);
        trisurf(S.(sn).mesh.triangles, X, Y, f([X, Y]));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
      end
    end
  end

  BBy = max(S.pos.numCellsY*S.pos.periodY, S.neg.numCellsY*S.neg.periodY);
  
  figure(1);
  xlim([-S.neg.numCellsX, S.pos.numCellsX]);
  ylim([-BBy, BBy]);

  figure(2);
  xlim([-S.neg.numCellsX, S.pos.numCellsX]);
  ylim([-BBy, BBy]);

  % figure(3);
  % x = [zeros(256, 1), linspace(-2, 2, 256)', zeros(256, 1)];
  % plot(x(:, 2), G(x));
  error;
end

%% Bilinear and linear forms
u = pdes.PDEObject; v = dual(u);
% gradu_gradv = @(muco) (muco * grad2(u)) * grad2(v);
% gradu_vecYv = @(muco) (muco * grad2(u)) * (EX * v);
% vecYu_gradv = @(muco) (muco * (EX * u)) * grad2(v);
% vecYu_vecYv = @(muco) (muco * (EX * u)) * (EX * v);
% u_v = @(rhoco) (rhoco * u) * v;

%% Compute FE elementary matrices
for Iside = 1:2
  sn = sides.name{Iside};
  dom = S.(sn).mesh.domain('volumic');

  S.(sn).mat_gradu_gradv = FEPack.pdes.Form.intg(dom, (S.(sn).mu * grad2(u)) * grad2(v));
  S.(sn).mat_gradu_vecYv = FEPack.pdes.Form.intg(dom, (S.(sn).mu * grad2(u)) * (EX * v));
  S.(sn).mat_vecYu_gradv = FEPack.pdes.Form.intg(dom, (S.(sn).mu * (EX * u)) * grad2(v));
  S.(sn).mat_vecYu_vecYv = FEPack.pdes.Form.intg(dom, (S.(sn).mu * (EX * u)) * (EX * v));
  S.(sn).mat_u_v = FEPack.pdes.Form.intg(dom, (S.(sn).rho * u) * v);
  S.(sn).mat_u_v_id = FEPack.pdes.Form.intg(dom, u * v);
end

%% Boundary conditions
for Iside = 1:2
  sn = sides.name{Iside};
  FourierIds = [0, floor(S.(sn).numNodesY/4)];
  
  S.(sn).BCstruct.spB0 = spaces.FourierBasis(S.(sn).mesh.domain('xmin'), FourierIds);
  S.(sn).BCstruct.spB1 = spaces.FourierBasis(S.(sn).mesh.domain('xmax'), FourierIds);

  S.(sn).BCstruct.BCu  = 1.0;
  S.(sn).BCstruct.BCdu = 0.0;
  S.(sn).BCstruct.representation = 'weak evaluation';
end 

%% Floquet-Bloch transform of auxiliary solutions
for Iside = 1:2
  sn = sides.name{Iside};
  ss = sides.sign{Iside};

  NFB = S.(sn).numFBpts;
  S.(sn).FBpts  = linspace(-pi/S.(sn).periodY, pi/S.(sn).periodY, NFB);
  S.(sn).sol    = cell(NFB, 1);
  S.(sn).DtN    = cell(NFB, 1);
  S.(sn).intsol = cell(NFB, 1);
  S.(sn).intDtN = cell(NFB, 1);
  S.(sn).intrhs = cell(NFB, 1);
  
  % FBT source term
  sourceFB = FEPack.tools.BlochTransform(S.(sn).mesh.points, S.(sn).FBpts', source, 2, S.(sn).periodY);

  % For the interior problem
  N = S.(sn).mesh.numPoints;
  spB0 = S.(sn).BCstruct.spB0;
  Nb = spB0.numBasis; 
  B0 = sparse(N, Nb);

  ecs = FEPack.pdes.EssentialConditions;           % BCs for interior problem
  ecs = ecs & ...
      ((u|S.(sn).mesh.domain('ymin')) -...
       (u|S.(sn).mesh.domain('ymax')) == 0.0) &... % Periodic BC
      ((u|S.(sn).mesh.domain('xmin')) == 0.0) &... % Dirichlet
      ((u|S.(sn).mesh.domain('xmax')) == 0.0);
  ecs.applyEcs;
  ecs.b = B0;

  % REMARQUE: Techniquement, le dernier (NFB) n'a 
  % pas besoin d'être calculé
  for idFB = 1:NFB
    fprintf('%d sur %d\n', idFB, NFB);
    FBvar = S.(sn).FBpts(idFB);

    % FE matrix
    % |||||||||
    AA = S.(sn).mat_gradu_gradv +...
         S.(sn).mat_vecYu_gradv * (+1i * FBvar) +...
         S.(sn).mat_gradu_vecYv * (-1i * FBvar) +...
         S.(sn).mat_vecYu_vecYv * (FBvar^2) +...
         S.(sn).mat_u_v * (-opts.omega^2);
    
    % Solve half-waveguide problem
    % ||||||||||||||||||||||||||||
    [~, S.(sn).sol{idFB}.E0,...
        S.(sn).sol{idFB}.E1,...
        S.(sn).sol{idFB}.R, ... 
        S.(sn).sol{idFB}.D, ...
        S.(sn).BCstruct,...
        S.(sn).DtN{idFB}] = PeriodicHalfGuideBVP(S.(sn).mesh, ss, 1, AA, S.(sn).BCstruct, S.(sn).numCellsX, opts);
    
    % Solve interior problems
    % |||||||||||||||||||||||
    TT = FEPack.pdes.Form.intg_TU_V(...
          S.(sn).mesh.domain('xmax'),...
          S.(sn).DtN{idFB},...
          S.(sn).BCstruct.spB0, 'projection');

    AAint = AA + TT;
    LLsource = S.(sn).mat_u_v_id * sourceFB(:, idFB);

    AA0int = ecs.P * AAint * ecs.P';
    LL0int = ecs.P * [LLsource, -AA * ecs.b];
    
    S.(sn).intsol{idFB} = [zeros(N, 1), ecs.b] + ecs.P' * (AA0int \ LL0int);
    S.(sn).intDtN{idFB} = spB0.invmassmat * B0' * AAint * S.(sn).intsol{idFB}(:, 2:end);
    S.(sn).intrhs{idFB} = spB0.invmassmat * B0' * AAint * S.(sn).intsol{idFB}(:, 1);
  end % for idFB
end % for Iside

clear AA AAint AA0int LLsource LL0int B0 ecs spB0 sourceFB N Nb;

%% Equation at the interface
% We use Lagrange P1 discretization
NFB  = S.pos.numFBpts;
NBZ  = S.pos.BCstruct.spB0.numBasis; 
NDOF = NBZ * (NFB - 1) + 1;
hpos = 2*pi / (S.pos.periodY * (NFB - 1));
S.pos.mat_DtN = sparse(NDOF, NDOF);

for Ielt = 1:NFB-1
  for Icell = 1:NBZ
    for Jcell = 1:NBZ
      mat_elem = [S.pos.DtN{Ielt}(Icell, Jcell) + 3*S.pos.DtN{Ielt+1}(Icell, Jcell),...
                  S.pos.DtN{Ielt}(Icell, Jcell) +   S.pos.DtN{Ielt+1}(Icell, Jcell);...
                  S.pos.DtN{Ielt}(Icell, Jcell) +   S.pos.DtN{Ielt+1}(Icell, Jcell),...
                3*S.pos.DtN{Ielt}(Icell, Jcell) +   S.pos.DtN{Ielt+1}(Icell, Jcell)];

      mat_elem = mat_elem * hpos/12;

      Iglob1 = Ielt + (Icell-1)*(NFB-1); 
      Iglob2 = Iglob1 + 1;
      Jglob1 = Ielt + (Jcell-1)*(NFB-1); 
      Jglob2 = Jglob1 + 1;

      if (Iglob1 >= 1 && Iglob1 <= NDOF)
        if (Jglob1 >= 1 && Jglob1 <= NDOF)
          S.pos.mat_DtN(Iglob1, Jglob1) = S.pos.mat_DtN(Iglob1, Jglob1) + mat_elem(1, 1);
        end

        if (Jglob2 >= 1 && Jglob2 <= NDOF)
          S.pos.mat_DtN(Iglob1, Jglob2) = S.pos.mat_DtN(Iglob1, Jglob2) + mat_elem(1, 2);
        end
      end

      if (Iglob2 >= 1 && Iglob2 <= NDOF)
        if (Jglob1 >= 1 && Jglob1 <= NDOF)
          S.pos.mat_DtN(Iglob2, Jglob1) = S.pos.mat_DtN(Iglob2, Jglob1) + mat_elem(2, 1);
        end

        if (Jglob2 >= 1 && Jglob2 <= NDOF)
          S.pos.mat_DtN(Iglob2, Jglob2) = S.pos.mat_DtN(Iglob2, Jglob2) + mat_elem(2, 2);
        end
      end
    end %
  end %
end


dofPts = -NBZ*pi/S.pos.periodY + (0:NDOF-1) * hpos;

NFBneg = S.neg.numFBpts;
NBZneg  = S.neg.BCstruct.spB0.numBasis; 
% NDOFneg = NBZneg * (NFBneg - 1) + 1;
% eltPts = sort([dofPts, -NBZneg*pi/S.neg.periodY + (0:NDOFneg-1) * hneg]);
% Nelt = numel(eltPts);



pYpos = S.pos.periodY;
pYneg = S.neg.periodY;

hneg = 2*pi / (pYneg * (NFBneg - 1));
S.neg.mat_DtN = sparse(NDOF, NDOF);

for Icell = 1:NBZneg
  for Jcell = 1:NBZneg
    for Ielt = 1:NFBneg-1
      Imod = Icell - 1 - (NBZneg-1)/2;
      Jmod = Jcell - 1 - (NBZneg-1)/2;

      Imin = floor((NBZ*pi/pYpos - (2*Imod+1)*pi/pYneg + (Ielt-1)*hneg)/hpos) + 1;
      Imax =  ceil((NBZ*pi/pYpos - (2*Imod+1)*pi/pYneg +     Ielt*hneg)/hpos) + 1;

      Jmin = floor((NBZ*pi/pYpos - (2*Jmod+1)*pi/pYneg + (Ielt-1)*hneg)/hpos) + 1;
      Jmax =  ceil((NBZ*pi/pYpos - (2*Jmod+1)*pi/pYneg +     Ielt*hneg)/hpos) + 1;
      
      disp([Imin, Imax, Jmin, Jmax])
      if (Imin <= NDOF && Imax >= 1 && Jmin <= NDOF && Jmax>= 1)
        Imin = max(Imin, 1); Imax = min(Imax, NDOF);
        Jmin = max(Jmin, 1); Jmax = min(Jmax, NDOF);

        % disp([Imin, Imax, Jmin, Jmax])
        pts = sort([-NBZ*pi/pYpos + ((Imin:Imax)-1)*hpos + 2*pi*Imod/pYneg,...
                    -NBZ*pi/pYpos + ((Jmin:Jmax)-1)*hpos + 2*pi*Jmod/pYneg]);
        DtN = S.neg.DtN{Ielt  }(Icell, Jcell) * (-pi/pYneg + Ielt*hneg - pts)/hneg +...
              S.neg.DtN{Ielt+1}(Icell, Jcell) * (pts + pi/pYneg - (Ielt-1)*hneg)/hneg;

        Ipts = floor((pts - 2*pi*Imod/pYneg + NBZ*pi/pYpos)/hpos) + 1;
        Jpts = floor((pts - 2*pi*Jmod/pYneg + NBZ*pi/pYpos)/hpos) + 1;

        WI = ((Imin:Imax)' == Ipts) .* (-NBZ*pi/pYpos + Ipts*hpos - pts + 2*pi*Imod/pYneg)/hpos +...
             ((Imin:Imax)' == Ipts+1) .* (pts - 2*pi*Imod/pYneg + NBZ*pi/pYpos - (Ipts-1)*hpos)/hpos;

        WJ = ((Jmin:Jmax)' == Jpts) .* (-NBZ*pi/pYpos + Jpts*hpos - pts + 2*pi*Jmod/pYneg)/hpos +...
             ((Jmin:Jmax)' == Jpts+1) .* (pts - 2*pi*Jmod/pYneg + NBZ*pi/pYpos - (Jpts-1)*hpos)/hpos;

        WIdiff  = WI(:, 2:end) - WI(:, 1:end-1);
        WI = WI(:, 1:end-1);
        
        WJdiff  = WJ(:, 2:end) - WJ(:, 1:end-1);
        WJ = WJ(:, 1:end-1);

        DtNdiff = diag((pts(2:end)-pts(1:end-1)) .* (DtN(2:end) - DtN(1:end-1)));
        DtN = diag(DtN(1:end-1));

        mat_elem = WI*DtN*WJ'...
                 + (WIdiff*DtN*WJ' + WI*DtNdiff*WJ' + WI*DtN*WJdiff')/2 ...
                 + (WIdiff*DtNdiff*WJ' + WIdiff*DtN*WJdiff' + WI*DtNdiff*WJdiff')/3 ...
                 + (WIdiff*DtNdiff*WJdiff')/4;

        S.neg.mat_DtN(Imin:Imax, Jmin:Jmax) = S.neg.mat_DtN(Imin:Imax, Jmin:Jmax) + mat_elem;
      end
    end %
  end %
end %




error;





error;

%% Positive half-space DtN
TFBphiPosVec = zeros(mesh2Dpos.domain('xmin').numPoints, numBasisInt, numFloquetPoints_pos);
for idI = 1:numBasisInt
  TFBphiPosVec(:, idI, :) = FEPack.tools.BlochTransform(mesh2Dpos.points(mesh2Dpos.domain('xmin').IdPoints, 2),...
                                           FloquetPoints_pos, @(x) int_basis_functions(x, idI),...
                                           1, period_pos, 1000);
end

% TFBphiPos = TFB_Lagrange_P1_1D(mesh2Dint, mesh2Dpos.points(mesh2Dpos.domain('xmin').IdPoints, 2), FloquetPoints_pos, period_pos);
TFBphiPos = cell(numFloquetPoints_pos, 1);
lambda_pos = zeros(numBasisInt);
wpos = (2*pi / period_pos) / (numFloquetPoints_pos - 1);

for idFB = 1:numFloquetPoints_pos
  TFBphiPos{idFB} = BCstruct_pos.spB0.FE_to_spectral * TFBphiPosVec(:, :, idFB);
  TFBlambdaPos{idFB} = BCstruct_pos.spB0.massmat * TFBlambdaPos{idFB};

  lambda_pos = lambda_pos + wpos * TFBphiPos{idFB}' * TFBlambdaPos{idFB} * TFBphiPos{idFB};
end

%% Negative half-space DtN
TFBphiNegVec = zeros(mesh2Dneg.domain('xmin').numPoints, numBasisInt, numFloquetPoints_neg);
for idI = 1:numBasisInt
  TFBphiNegVec(:, idI, :) = FEPack.tools.BlochTransform(mesh2Dneg.points(mesh2Dneg.domain('xmin').IdPoints, 2),...
                                           FloquetPoints_neg, @(x) int_basis_functions(x, idI),...
                                           1, period_neg, 1000);
end

% TFBphiNeg = TFB_Lagrange_P1_1D(mesh2Dint, mesh2Dneg.points(mesh2Dneg.domain('xmin').IdPoints, 2), FloquetPoints_neg, period_neg); 
TFBphiNeg = cell(numFloquetPoints_neg, 1);
lambda_neg = zeros(numBasisInt);
wneg = (2*pi / period_neg) / (numFloquetPoints_neg - 1);

for idFB = 1:numFloquetPoints_neg
  TFBphiNeg{idFB} = BCstruct_neg.spB0.FE_to_spectral * TFBphiNegVec(:, :, idFB);
  TFBlambdaNeg{idFB} = BCstruct_neg.spB0.massmat * TFBlambdaNeg{idFB};

  lambda_neg = lambda_neg + wneg * TFBphiNeg{idFB}' * TFBlambdaNeg{idFB} * TFBphiNeg{idFB};
end

%% Solve the integral equation on the interface
mat_G_v_int = spBint.projmat * Gint(mesh2Dint.points);
% -------------------------------------------------------------------- %
% The minus sign comes from the definition of the Lambda when they are %
% the DtN operators (see PeriodicHalfGuideBVP.m)                       %
% -------------------------------------------------------------------- %
trace_solution = -(lambda_pos + lambda_neg) \ mat_G_v_int;             %
% -------------------------------------------------------------------- %

%% Deduce the FB transform of the solution
% FB transform of the solution's trace with period_pos
TFB_solution_pos = cell(numFloquetPoints_pos, 1);
for idFB = 1:numFloquetPoints_pos

  TFB_solution_pos{idFB} = zeros(mesh2Dpos.numPoints, numCellsSemiInfinite_pos);
  R0Phi = TFBphiPos{idFB} * trace_solution;
  R1Phi = sol_pos_data{idFB}.D * R0Phi;

  for idCell = 0:numCellsSemiInfinite_pos-1
    % Compute the solution in the current cell
    TFB_solution_pos{idFB}(:, idCell + 1) = sol_pos_data{idFB}.E0 * R0Phi + sol_pos_data{idFB}.E1 * R1Phi;

    % Update
    R0Phi = sol_pos_data{idFB}.R * R0Phi;
    R1Phi = sol_pos_data{idFB}.D * R0Phi;
  end

end

% FB transform of the solution's trace with period_neg
TFB_solution_neg = cell(numFloquetPoints_neg, 1);
for idFB = 1:numFloquetPoints_neg

  TFB_solution_neg{idFB} = zeros(mesh2Dneg.numPoints, numCellsSemiInfinite_neg);
  R0Phi = TFBphiNeg{idFB} * trace_solution;
  R1Phi = sol_neg_data{idFB}.D * R0Phi;

  for idCell = 0:numCellsSemiInfinite_neg-1
    % Compute the solution in the current cell
    TFB_solution_neg{idFB}(:, idCell + 1) = sol_neg_data{idFB}.E0 * R0Phi + sol_neg_data{idFB}.E1 * R1Phi;

    % Update
    R0Phi = sol_neg_data{idFB}.R * R0Phi;
    R1Phi = sol_neg_data{idFB}.D * R0Phi;
  end

end

%% Compute the half-space solution cell by cell
% Positive side
numCells = [1 1];
numCells(infiniteDirection) = 2*numCellsInfinite_pos;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_pos;
Nu = prod(numCells);
[I1, I2] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2]; % 2-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite_pos' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period_pos) ./ (numFloquetPoints_pos - 1)); % 1-by-1
U.positive = zeros(mesh2Dpos.numPoints, Nu);

for idFB = 1:numFloquetPoints_pos
  FloquetVar = FloquetPoints_pos(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh2Dpos.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period_pos); % 1-by-Nu
  U_TFB = TFB_solution_pos{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U.positive = U.positive + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U.positive = U.positive * sqrt(period_pos / (2*pi));

% Negative side
numCells = [1 1];
numCells(infiniteDirection) = 2*numCellsInfinite_neg;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_neg;
Nu = prod(numCells);
[I1, I2] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2]; % 2-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite_neg' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period_neg) ./ (numFloquetPoints_neg - 1)); % 1-by-1
U.negative = zeros(mesh2Dneg.numPoints, Nu);

for idFB = 1:numFloquetPoints_neg
  FloquetVar = FloquetPoints_neg(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh2Dneg.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period_neg); % 1-by-Nu
  U_TFB = TFB_solution_neg{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U.negative = U.negative + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U.negative = U.negative * sqrt(period_neg / (2*pi));

%% Plot the solution
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
for idS = 1:numCellsSemiInfinite_pos
  for idI = 1:2*numCellsInfinite_pos
    Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite_pos], idS, idI);
    X = mesh2Dpos.points(:, 1) + (idS - 1);
    Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite_pos - 1) * period_pos;

    trisurf(mesh2Dpos.triangles, X, Y, real(U.positive(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end
%
for idS = 1:numCellsSemiInfinite_neg
  for idI = 1:2*numCellsInfinite_neg
    Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite_neg], idS, idI);
    X = mesh2Dneg.points(:, 1) - (idS - 1);
    Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite_neg - 1) * period_neg;
    trisurf(mesh2Dneg.triangles, X, Y, real(U.negative(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

axis([-4, 4, -4, 4]);
% xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);

