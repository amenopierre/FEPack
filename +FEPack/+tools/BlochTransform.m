function TFBfun = BlochTransform(x, k, fun, directions, Lcell, Ntrunc, BlochType)
  % function TFBfun = BlochTransform(x, k, fun)
  % returns the Floquet-Bloch transform of a function along a given direction.
  % More precisely,
  %
  %  TFBfun(x) = prod(Lcell_j)/(2*pi)^(d) *
  %                         sum_{n in Z^d} fun(x + nI) exp(-i*<k, xI + n>),
  %
  % where nI_i = ni and xI_i = xi if i is a FB direction and 0 otherwise.
  %
  % INPUTS: * x, a Nx-by-3 sized vector containing the points.
  %         * k, a Nk-by-d vector which is Floquet variable.
  %         * fun, the function. fun should take a N-by-3 matrix in argument,
  %           and return a N-sized vector.
  %         * directions, a d-sized vector containing the directions along
  %           which the Floquet Bloch transform should be applied.
  %         * Lcell, a d-sized vector containing the sizes of the cells
  %         * Ntrunc, a d-sized vector containing the truncation sizes
  %           for the sum.
  %         * BlochType, a string indicating the Floquet-Bloch transform,
  %           either 'periodic' or 'quasiperiodic'.
  %
  % OUTPUTS: * TFBfun, a Nx-by-Nk matrix

  directions = unique(directions(:))';
  dFB = length(directions);
  x = [x, zeros(size(x, 1), 3-size(x, 2))];

  % Default variables and preliminary verifications
  % ///////////////////////////////////////////////
  if (nargin < 5)
    Lcell = ones(1, dFB);
  end

  if (nargin < 6)
    Ntrunc = 100*ones(1, dFB);
  end
  
  if (nargin < 7)
    BlochType = 'periodic';
  end

  if (length(Ntrunc) ~= length(directions))
    % Ntrunc and directions should have the same length
    error('Les vecteurs Ntrunc et directions doivent être de même longueur.');
  end

  if (length(Lcell) ~= length(directions))
    % Lcell and directions should have the same length
    error('Les vecteurs Lcell et directions doivent être de même longueur.');
  end

  if (size(k, 2) ~= length(directions))
    % k should have as much columns as directions
    error('k doit avoir autant de colonnes que le nombre de directions.')
  end

  for idI = 1:dFB
    if min((1:3) ~= directions(idI))
      % The components along which the FB transform is applied are >= 1 and <= 3
      error('La TFB ne peut pas être appliquée dans la direction %d.', directions(idI));
    end
  end
  
  if ~(strcmpi(BlochType, 'periodic') || strcmpi(BlochType, 'quasiperiodic'))
    % The type of Floquet-Bloch transform should be either 'periodic' or
    % 'quasiperiodic'
    error(['Les valeurs possibles pour BlochType sont ''periodic'' et ',...
           '''quasiperiodic''. L''argument ''%s'' n''est pas possible.'], BlochType);
  end

  % Construct the set of indices
  % ////////////////////////////
  sumVars = cell(dFB, 1);
  for idI = 1:dFB
    sumVars{idI} = -Ntrunc(idI):Ntrunc(idI);
  end

  Nsum = prod(2*Ntrunc + 1);
  Dsum = [1 1 1];
  Dsum(directions) = 2*Ntrunc+1;
  [I1, I2, I3] = ind2sub(Dsum, 1:Nsum);
  sumIds = [I1; I2; I3];
  sumIds = sumIds(directions, :);

  NLcell = zeros(Nsum, 3);
  for idI = 1:dFB
    NLcell(:, directions(idI)) = Lcell(idI) * sumVars{idI}(sumIds(idI, :))';
  end

  % Compute the periodic Floquet Bloch transform
  % ////////////////////////////////////////////
  Nx = size(x, 1);

  x_plus_NL = kron(ones(Nsum, 1), x) + kron(NLcell, ones(Nx, 1)); % Nsum*Nx-by-3
  F_x_plus_NL = reshape(fun(x_plus_NL), Nx, Nsum);                % Nx-by-Nsum

  k_dot_NL = NLcell(:, directions) * k.';  % Nsum-by-Nk
  exp_k_dot_NL = exp(-1i * k_dot_NL);      % Nsum-by-Nk

  k_dot_x = x(:, directions) * k.';        % Nx-by-Nk
  exp_k_dot_x = exp(-1i * k_dot_x);        % Nx-by-Nk
  
  TFBfun = prod(sqrt(Lcell/(2*pi))) * exp_k_dot_x .* (F_x_plus_NL * exp_k_dot_NL);
  
  % Switch to the quasiperiodic Floquet-Bloch transform if needed
  % /////////////////////////////////////////////////////////////
  if strcmpi(BlochType, 'quasiperiodic')
    
    TFBfun = TFBfun .* exp(1i * k_dot_x);
    
  end
  
end
