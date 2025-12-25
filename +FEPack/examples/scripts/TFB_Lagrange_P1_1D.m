function TFBphi = TFB_Lagrange_P1_1D(mesh, x, k, period, BlochType)

  if (nargin < 5)
    BlochType = 'periodic';
  end
  x = x(:);
  k = k(:);
  Nx = size(x, 1);
  Nint = mesh.numPoints - 2;
  Nk = length(k);
  TFBvec = zeros(Nx, Nint, Nk);

  for idB =  1:Nint

    nodes = mesh.points(idB + 1, 1);
    nodes_plus_one = mesh.points(idB + 2, 1);
    nodes_minus_one = mesh.points(idB, 1);

    phi = @(P) ((P(:, 1) - nodes_minus_one) ./ (nodes - nodes_minus_one)) .* (P(:, 1) <= nodes & P(:, 1) > nodes_minus_one) + ...
               ((nodes_plus_one -  P(:, 1)) ./ (nodes_plus_one -  nodes)) .* (P(:, 1) <= nodes_plus_one &  P(:, 1) > nodes);
    %
    TFBvec(:, idB, :) = FEPack.tools.BlochTransform(x, k, phi, 1, period, 1000, BlochType);

  end
  
  TFBphi = cell(Nk, 1);
  for idFB = 1:Nk
    TFBphi{idFB} = TFBvec(:, :, idFB);
  end


  
  
  
  
  
  
  
  
  
  
  % if (nargin < 5)
  %   BlochType = 'periodic';
  % end
  
  % Nx = size(x, 1);
  % Nint = mesh.numPoints - 2;

  % x_kron = x * ones(1, Nint);
  % int_nodes_kron = ones(Nx, 1) * mesh.points(2:end-1, 1).';
  % int_nodes_kron_plus_one = ones(Nx, 1) * mesh.points(3:end, 1).';
  % int_nodes_kron_minus_one = ones(Nx, 1) * mesh.points(1:end-2, 1).';
  
  % x_kron = x_kron(:);
  % int_nodes_kron = int_nodes_kron(:);
  % int_nodes_kron_plus_one = int_nodes_kron_plus_one(:);
  % int_nodes_kron_minus_one = int_nodes_kron_minus_one(:);

  % % Indices for which phi(x + n*p) is non-null
  % idMIN = ceil((int_nodes_kron_minus_one - x_kron) / period);
  % idMAX = floor((int_nodes_kron_plus_one - x_kron) / period);
  % numTerms = max(max(idMAX - idMIN + 1, 0));

  % x_plus_NL = (x_kron + idMIN * period) * ones(1, numTerms) + ones(Nx * Nint, 1) * ((0:numTerms-1) * period);

  % phi_x_plus_NL = zeros(Nx * Nint, numTerms);
  % idA = (x_kron >= int_nodes_kron) & (x_kron < int_nodes_kron_plus_one);
  % idB = (x_kron >= int_nodes_kron_minus_one) & (x_kron < int_nodes_kron);

  % phi_x_plus_NL(idA) = (int_nodes_kron_plus_one(idA) - x_plus_NL(idA)) ./ (int_nodes_kron_plus_one(idA) - int_nodes_kron(idA)); 
  % phi_x_plus_NL(idB) = (x_plus_NL(idB) - int_nodes_kron_minus_one(idB)) ./ (int_nodes_kron(idB) - int_nodes_kron_minus_one(idB)); 

  % %% Compute the FB transform
  % Nk = length(k);
  % TFBphi = cell(Nk, 1);

  % for idK = 1:Nk

  %   exp_k_dot_NL = exp(-1i * k(idK) * x_plus_NL);

  %   TFBphi{idK} = sqrt(period / (2 * pi)) * sum(phi_x_plus_NL .* exp_k_dot_NL, 2);
  %   TFBphi{idK} = reshape(TFBphi{idK}, Nx, Nint);
    
  % end

  % if strcmpi(BlochType, 'quasiperiodic')
  %    for idK = 1:Nk

  %     conj_exp_k_dot_x = exp(1i * k(idK) * x_kron);

  %     TFBphi{idK} = reshape(TFBphi{idK}(:) .* conj_exp_k_dot_x, Nx, Nint);

  %   end    
  % end


end
