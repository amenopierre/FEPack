%> @file Form.m
%> @brief Contains the pdes.Form class.
% =========================================================================== %
%> @brief class for PDE object
% =========================================================================== %
classdef Form < FEPack.FEPackObject
  % FEPack.pdes.Form < FEPack.FEPackObject

  properties (SetAccess = public)

    %> @brief Indicates if linear or bilinear form
    is_bilinear = 1;

    %> @brief Coefficients associated to primal variable
    alpha_u = [];

    %> @brief Coefficients associated to dual variable
    alpha_v = [];

    %> @brief Multiplicative coefficient
    fun = @(x) ones(size(x, 1), 1);

  end

  methods (Static)

    % Shape functions
    function phis = shapeFunctions(P, dimension, alpha)
      % SHAPEFUNCTIONS Expression of derivatives of Lagrange P1 shape functions.
      % phis = (P, d, alpha) where P is a point and alpha a (d+1)-vector,
      % returns the linear combination
      %
      %   alpha_1 * phi + \sum_{i = 1}^d alpha_{i+1} * (d_xi phi)
      %
      % where phi is a Lagrange P1 shape function.
      %
      % INPUTS: * P, the points in which the shape functions are evaluated.
      %           P is of size N-by-d
      %         * dimension, the dimension of the points (1, 2, or 3)
      %         * alpha, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the derivatives of the shape functions.
      %
      % OUTPUTS: * phis, (d+1)-by-N matrix containing the shape functions

      d = dimension;

      N = size(P, 1);
      phis = zeros(N, d+1);

      phis(:, 1) = alpha(1) * (1 - P * ones(d, 1)) - ones(N, d) * alpha(2:d+1)';
      phis(:, 2:end) = alpha(1) * P + kron(ones(N, 1), alpha(2:d+1));

    end

    % Elementary matrices
    function Aelem = mat_elem(P, domDim, volDim, alpha_u, alpha_v, fun, quadRule)
      % Compute elementary matrix whose components are given by
      %
      %  \int_{T} fun(x) * [alpha_1 * u + \sum_{i = 1}^d alpha_{i+1} * (d_xi u)]
      %                  * [ beta_1 * v + \sum_{i = 1}^d  beta_{i+1} * (d_xi v)]
      %
      % where T is a given element, and where u and v are the Lagrange P1 basis
      % functions associated to the nodes of the element T.
      %
      % INPUTS: * P, 3x(d+1) matrix containing the coordinates of the element's
      %           vertices, where d is the domain dimension.
      %         * domDim, the dimension of the domain.
      %           WARNING: not to be confused with the dimension of the
      %                    volumic geometry. For instance, for a segment,
      %                    domDim = 1.
      %         * volDim, the dimension of the volumic geometry
      %         * alpha_u, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the unknown and its derivatives.
      %         * alpha_v, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the test function and its derivatives.
      %         * fun (function handle) the coefficient in the integral
      %           WARNING: fun must take in argument a nx3 matrix, and return
      %                    a nx1 vector.
      %         * quadRule (QuadratureObject) is optional. For more information,
      %           see +tools/QuadratureObject.m
      %
      % OUTPUTS: Aelem, a (domDim+1)-by-(domDim+1) matrix.

      % Each point must have volDim coordinates at least and 3 coordinates at most
      if ((size(P, 1) < volDim) || (size(P, 1) > 3))
        error(['Les points doivent un nombre de coordonnées égal à 3 au ',...
               'plus, et à la dimension volumique (', int2str(volDim), ') au moins.']);
      end

      % The domain dimension should match the number of points in the element
      if (domDim + 1 ~= size(P, 2))
        error(['La dimension du domaine + 1 (', int2str(domDim + 1), ') et ',...
               'le nombre de points (', int2str(size(P, 2)), ') de l''élément ',...
               'doivent coincider.']);
      end

      % The domain dimension should not exceed the volumic dimension
      if (domDim > volDim)
        error(['La dimension du domaine (', int2str(domDim), ') ne peut ',...
               'pas être supérieure à la dimension volumique (', int2str(volDim), ')']);
      end

      % Preliminary adjustments and default values
      Lu = length(alpha_u); alpha_u = [alpha_u(:); zeros(4-Lu, 1)]; % Fill alpha_u with zeros
      Lv = length(alpha_v); alpha_v = [alpha_v(:); zeros(4-Lv, 1)]; % Fill alpha_v with zeros
      dP = size(P); P = [P; zeros(3-dP(1), dP(2))];

      if (nargin < 6)
        fun = @(x) ones(size(x, 1), 1);
      end
      if (nargin < 7)
        quadRule = FEPack.tools.QuadratureObject(domDim);
      end
      Xquad = quadRule.points;
      Wquad = quadRule.weights;

      % Map to reference element
      mapToRel.A = -P(1:volDim, 1) + P(1:volDim, 2:end);
      mapToRel.B =  P(1:volDim, 1);
      switch (domDim)
      case 1
        mapToRel.J = sqrt((P(:, 2) - P(:, 1))' * (P(:, 2) - P(:, 1)));
      case 2
        Tu = cross(P(:, 2) - P(:, 1), P(:, 3) - P(:, 1));
        mapToRel.J = sqrt(Tu' * Tu);
      case 3
        mapToRel.J = abs(det(mapToRel.A));
      end

      % Weighted function
      weightedFun = Wquad .* fun((mapToRel.A * Xquad + mapToRel.B).').';

      % Shape functions
      coeffs.u = alpha_u;
      coeffs.v = alpha_v;
      fieldnames = ['u'; 'v'];

      for idI = 1:2
        beta = coeffs.(fieldnames(idI));

        if (domDim == volDim)
          % 0-order term and partial derivatives
          alpha = [beta(1), (mapToRel.A \ beta(2:domDim+1)).'];
        else
          % 0-order term and tangential derivatives
          % WARNING: this should be used with caution (at least validation
          % needed)
          alpha = [beta(1), beta(2:domDim+1).' ./ sqrt(diag(mapToRel.A' * mapToRel.A)).' ];
        end

        phis.(fieldnames(idI)) = FEPack.pdes.Form.shapeFunctions(Xquad.', domDim, alpha);
      end

      % Compute the elementary matrix
      Aelem = mapToRel.J * (phis.v' * diag(weightedFun) * phis.u);

    end

    % Matrix assembly
    function Aglob = assembleFEmatrices(domain, Aloc)
      % ASSEMBLEFEMATRICES Finite elements matrix assembly
      % Aglob = ASSEMBLEFEMATRICES(domain, Aloc) where domain is a
      % FEPack.meshes.FEDomain object, and where Aloc is a function handle
      % computes a sparse matrix that results from the assembly of local
      % matrices computed on the elements of domain, via Aloc.
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object, the domain for which
      %           the assembly process is performed.
      %         * Aloc a function handle. Given a nx3 matrix corresponding
      %           containg the coordinates of the n points of an element,
      %           Aloc should return the n-by-n associated elementary matrix.
      %
      % OUTPUTS: * Aglob, a N-by-N matrix, where N is the number of DOFs.

      N = domain.mesh.numPoints;   % Number of degrees of freedom
      domDim = domain.dimension;
      dd = (domDim + 1) * (domDim + 1);

      II = zeros(domain.numElts*dd, 1);
      JJ = zeros(domain.numElts*dd, 1);
      VV = zeros(domain.numElts*dd, 1);

      index_II = kron(ones(domDim + 1, 1), (1:domDim + 1).');
      index_JJ = kron((1:domDim + 1).', ones(domDim + 1, 1));

      for ielts = 1:domain.numElts

        % Nodes composing to the element
        P = domain.mesh.points(domain.elements(ielts, :), :).';

        % Elementary matrix associated to the element
        Aelem = Aloc(P);

        % Save the elementary matrix
        index = (dd*(ielts-1)+1):(dd*(ielts-1)+dd);
        II(index) = domain.elements(ielts, index_II);
        JJ(index) = domain.elements(ielts, index_JJ);
        VV(index) = Aelem(:);

      end

      Aglob = sparse(II, JJ, VV, N, N);

    end

    % Global matrices
    function AA = global_matrix(varargin)

      % function AA = GLOBAL_MATRIX(domain, alpha_u, alpha_v, fun, quadRule).
      %
      % If U = [u, (d_x1 u),..., (d_xN u)] and V = [v, (d_x1 v),..., (d_xN v)],
      % computes the FE matrix associated to the bilinear form
      %
      %  \int_{domain} [fun(x) * (Au * U)] * (Av * V),
      %
      % where
      %   * U = [u, (d_x1 u),..., (d_xN u)], V = [v, (d_x1 v),..., (d_xN v)],
      %   * Au and Av are p-by-(N+1) matrices.
      %
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object, the domain on which
      %           the integrals are evaluated.
      %         * alpha_u, a p-by-(d+1)-vector that contains the coefficients
      %           of the linear combination of the unknown and its derivatives.
      %         * alpha_v, a p-by-(d+1)-vector that contains the coefficients
      %           of the linear combination of the test function and its derivatives.
      %         * fun (function handle) the coefficient in the integral
      %           WARNING: fun must take in argument a n-by-3 matrix, and return
      %                    a n-by-p matrix. An example of such function would be
      %
      %                 matfun = @(P) [cos(P(:, 1)), sin(P(:, 2));...
      %                                exp(P(:, 1)), zeros(size(P, 1), 1)]
      %
      %         * quadRule (QuadratureObject) is optional. For more information,
      %           see +tools/QuadratureObject.m
      %
      % OUTPUTS: * AA, a N-by-N matrix, where N is the number of DOFs.

      dom = varargin{1};
      alpha_u = varargin{2};
      alpha_v = varargin{3};
      Ncoo = size(alpha_u, 1);

      % Computing tangential derivative on surfacic domain is not allowed (yet)
      if (((dom.dimension < dom.mesh.dimension) || (dom.dimension == 0)) && ...
           (~isempty(find(alpha_u(2:end) == 1, 1)) || ~isempty(find(alpha_v(2:end) == 1, 1))))
        error(['Calculer l''intégrale d''une dérivée tangentielle sur un ',...
               'domaine surfacique n''est pas autorisé. Pour lever cette ',...
               'erreur, revoir le calcul des matrices élémentaires.']);
      end

      % Make sure alpha_u and alpha_v have compatible size
      % More precisely, fun(P)*alpha_u and alpha_v must have the same size.
      if (((nargin < 4) &&  min(size(alpha_u) ~= size(alpha_v))) ||...
          ((nargin > 3) &&    ((size(varargin{4}([0, 0, 0]), 2) ~= size(alpha_u, 1)) ||...
                               (size(varargin{4}([0, 0, 0]), 1) ~= size(alpha_v, 1)))))
        error(['Les tailles de alpha_u, alpha_v, et de la sortie de fun ',...
               'doivent être de telle sorte que fun(P)*alpha_u soit de la ',...
               'même taille que alpha_v']);
      end

      AA = sparse(0);

      if (nargin < 4)

        % fun = 1 by default
        for Icoo = 1:Ncoo
          Aloc = @(P) FEPack.pdes.Form.mat_elem(P, dom.dimension, dom.mesh.dimension, alpha_u(Icoo, :), alpha_v(Icoo, :));

          AA = AA + FEPack.pdes.Form.assembleFEmatrices(dom, Aloc);
        end

      else

        for Icoo = 1:Ncoo
          for Jcoo = 1:Ncoo
            % Extract the desired component of the matrix function
            EI = @(P) sparse(1:size(P,1), ((Icoo-1)*size(P,1)+1):(Icoo*size(P,1)), 1, size(P,1), Ncoo*size(P,1));
            EJ = @(P) sparse(Jcoo, 1, 1, Ncoo, 1);

            % Elementary matrix
            fun = @(P) EI(P) * varargin{4}(P) * EJ(P); % (Icoo, Jcoo)-component of matrix function
            Aloc = @(P) FEPack.pdes.Form.mat_elem(P, dom.dimension, dom.mesh.dimension, alpha_u(Icoo, :), alpha_v(Icoo, :), fun, varargin{5:end});

            % Update the matrix
            AA = AA + FEPack.pdes.Form.assembleFEmatrices(dom, Aloc);
          end
        end

      end

    end

    % More user-friendly aliases
    function AA = intg_U_V(varargin)
      % function AA = intg_U_V(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, 1, 1, varargin{2:end});
    end

    function AA = intg_U_DxV(varargin)
      % function AA = intg_U_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [1 0], [0 1], varargin{2:end});
    end

    function AA = intg_U_DyV(varargin)
      % function AA = intg_U_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [1 0 0], [0 0 1], varargin{2:end});
    end

    function AA = intg_U_DzV(varargin)
      % function AA = intg_U_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [1 0 0 0], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_DxU_V(varargin)
      % function AA = intg_DxU_V(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 1], [1 0], varargin{2:end});
    end

    function AA = intg_DxU_DxV(varargin)
      % function AA = intg_DxU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 1], [0 1], varargin{2:end});
    end

    function AA = intg_DxU_DyV(varargin)
      % function AA = intg_DxU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 1 0], [0 0 1], varargin{2:end});
    end

    function AA = intg_DxU_DzV(varargin)
      % function AA = intg_DxU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 1 0 0], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_DyU_V(varargin)
      % function AA = intg_DyU_V(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 1], [1 0 0], varargin{2:end});
    end

    function AA = intg_DyU_DxV(varargin)
      % function AA = intg_DyU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 1], [0 1 0], varargin{2:end});
    end

    function AA = intg_DyU_DyV(varargin)
      % function AA = intg_DyU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 1], [0 0 1], varargin{2:end});
    end

    function AA = intg_DyU_DzV(varargin)
      % function AA = intg_DyU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 1 0], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_DzU_V(varargin)
      % function AA = intg_DzU_V(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 0 1], [1 0 0 0], varargin{2:end});
    end

    function AA = intg_DzU_DxV(varargin)
      % function AA = intg_DzU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 0 1], [0 1 0 0], varargin{2:end});
    end

    function AA = intg_DzU_DyV(varargin)
      % function AA = intg_DzU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 0 1], [0 0 1 0], varargin{2:end});
    end

    function AA = intg_DzU_DzV(varargin)
      % function AA = intg_DzU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, [0 0 0 1], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_gradU_gradV(varargin)
      % function AA = intg_gradU_gradV(domain, matfun, quadRule)
      d = varargin{1}.dimension;
      alpha = [zeros(d, 1), eye(d)];
      AA = FEPack.pdes.Form.global_matrix(varargin{1}, alpha, alpha, varargin{2:end});
    end

    function AA = intg_TU_V(domain, Tmat, spectralB, representation)
      % function AA = INTG_TU_V(domain, Tmat, space, representation)
      % Computes the FE matrix associated to the integral
      %
      %       \intg_domain (Tu) * v
      %
      % where T is an operator applied to u. T is represented by a matrix by
      % means of a spectral basis (phi_k)_k. The components of the matrix are
      % the Tmat_{k,l} which are either given by
      %
      %       (1) weak evaluation: Tmat_{k,l} = <Tphi_l, phi_k>; or
      %
      %       (2) projection and representation:
      %             Proj(Tphi_l) = Tmat_{1,l} phi_1 + .... + Tmat_{N, l} phi_N
      %          where Proj(Tphi_l) is the projection of phi_l in the (phi_k)_k.
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object, the domain on which
      %           the integrals are evaluated.
      %         * Tmat, a matrix that represents the operator applied to the
      %           unknown.
      %         * spectralB, SpectralBasis object.
      %         * representation, a string between 'weak evaluation' and
      %           'projection', which specifies the definition of T.
      %
      % OUTPUTS: * AA, a N-by-N matrix, where N is the number of DOFs.

      if isa(Tmat, 'function_handle')

        % The operator is a multiplication by a function
        AA = FEPack.pdes.Form.intg_U_V(domain, Tmat);

      elseif (length(Tmat) == 1)

        % Trivial case of multiplication by a scalar
        AA = Tmat * FEPack.pdes.Form.intg_U_V(domain);
        
      elseif ~(size(Tmat) == spectralB.numBasis)

        error('Si T est une matrice, alors sa taille doit être égale au nombre de fonctions de base spectrale.');

      else

        if strcmpi(representation, 'projection')

          % Deduce the weakly evaluated matrix from the projected one
          Tmat = spectralB.massmat * Tmat;

        elseif ~strcmpi(representation, 'weak evaluation')

          % Only 'projection' and 'weak evaluation' are allowed
          error(['La variable evaluation ne peut valoir que ''weak ',...
                 'evaluation'' ou ''projection''.']);

        end

        % If not done already, compute the matrices associated to
        % the spectral basis
        if isempty(spectralB.massmat)
          spectralB.computeBasisMatrices(0);
        end

        % Deduce the matrix
        Proj = spectralB.massmatInv * spectralB.projmat.';

        N = domain.mesh.numPoints;
        AA = sparse(N, N);
        AA(domain.IdPoints, domain.IdPoints) = Proj' * Tmat * Proj;

      end

    end

    function AA = intg(varargin)
      % function AA = intg(domain, aLF, quadRule)
      %
      % Compute a matrix or a vector associated to a bilinear or linear form.

      aLF = varargin{2};

      if (isa(aLF, 'FEPack.pdes.LinOperator') && aLF.is_dual)
        % Linear form
        alpha_u = [1 0 0 0];
        alpha_v = aLF.alpha;
      elseif isa(aLF, 'FEPack.pdes.Form')
        % Bilinear form
        alpha_u = aLF.alpha_u;
        alpha_v = aLF.alpha_v;
      else
        error(['La fonction intg ne peut être appliquée qu''à un produit ',...
               'd''opérateurs (LinOperator) sur inconnue et sur fonction test ',...
               '(forme bilinéaire), ou à un opérateur sur fonction test (forme ',...
               'linéaire)']);
      end

      AA = FEPack.pdes.Form.global_matrix(varargin{1}, alpha_u, alpha_v, aLF.fun, varargin{3:end});

      % For linear forms, deduce the vector
      if (isa(aLF, 'FEPack.pdes.LinOperator') && aLF.is_dual)
        AA = AA * ones(size(AA, 2), 1);
      end

    end

  end

  methods

    function aLFres = plus(aLFA, aLFB)
      aLFres = copy(aLFA);
      aLFres.alpha_u = [aLFA.alpha_u; aLFB.alpha_u];
      aLFres.alpha_v = [aLFA.alpha_v; aLFB.alpha_v];
      aLFres.fun = @(P) blkdiag(aLFA.fun(P), aLFB.fun(P));
    end

    function aLFres = mtimes(T, aLF)
      aLFres = copy(aLF);
      aLFres.fun = @(P) T * aLF.fun(P);
    end

    function aLFres = uminus(aLF)
      aLFres = (-1) * aLF;
    end

    function aLFres = minus(aLFA, aLFB)
      aLFres = aLFA + (-1) * aLFB;
    end

  end

end
