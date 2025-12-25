%> @file HoneycombObject.m
%> @brief Contains the +applications.HoneycombObject class.
% =========================================================================== %
%> @brief class for examples of honeycomb potentials
% =========================================================================== %
classdef HoneycombObject < FEPack.FEPackObject
  % FEPack.applications.HoneycombObject < FEPack.FEPackObject

  properties (SetAccess = public)

    % Expression of the honeycomb lattice potentials
    V = [];
    W = [];
    A = [];
    
    % Center of lattice
    center = [];

    % Lattice vectors
    vecPer1 = [];
    vecPer2 = [];

    % Dual lattice vectors
    dualVec1 = [];
    dualVec2 = [];

    % High symmetry quasi-momenta
    highSymK = [];
    brillouinVerts = [];

    % Edge
    edge = [];

    % Domain wall
    kappa = [];

  end

  methods

    function obj = HoneycombObject(typeV, typeW, typeA, funV, funW, funA)

      % is_even (boolean): is the potential even?
      % is_real (boolean): is the potential real-valued?
      % type (string): options are 'atomic', 'optical', 'trigonometric', 'custom', and 'none'
      % fun (function handle): only need if the type is 'custom'; gives a user-specified expression of the potential

      % Center of lattice
      obj.center = [0; 0];

      % Lattice vectors
      obj.vecPer1 = [0.5*sqrt(3);  0.5];
      obj.vecPer2 = [0.5*sqrt(3); -0.5];

      % Dual lattice vectors
      obj.dualVec1 = (4*pi/sqrt(3)) * [0.5;  0.5*sqrt(3)];
      obj.dualVec2 = (4*pi/sqrt(3)) * [0.5; -0.5*sqrt(3)];

      % High symmetry quasi-momenta
      obj.highSymK = (obj.dualVec1 - obj.dualVec2) / 3; % High-symmetry quasi-momentum
      obj.brillouinVerts = [obj.highSymK,...
                           -obj.highSymK + obj.dualVec1,...
                            obj.highSymK + obj.dualVec2,...
                           -obj.highSymK,...
                            obj.highSymK - obj.dualVec1,...
                           -obj.highSymK - obj.dualVec2];
      
      % V
      if strcmpi(typeV, 'atomic')

        obj.V = @(x) FEPack.tools.atomicPotential(x, obj.vecPer1, obj.vecPer2);

      elseif strcmpi(typeV, 'optical')

        obj.V = @(x) FEPack.tools.opticalPotential(x, obj.dualVec1, obj.dualVec2, true);

      elseif strcmpi(typeV, 'trigonometric')

        obj.V = @(x) obj.trigonometricPolynomial(x, true);

      elseif strcmpi(typeV, 'custom')

        if (nargin < 3)
          error(['Option', typeV, ' was selected: a function handle has to be provided.']);
        end

        obj.V = funV;

      elseif strcmpi(typeV, 'none')

        obj.V = [];
      
      else

        error(['Option ', typeV, ' unrecognized; only options are ''atomic'', ''optical'', ''trigonometric'', ''custom''']);

      end

      % A
      if strcmpi(typeA, 'atomic')

        obj.A = @(x) FEPack.tools.atomicPotential(x, obj.vecPer1, obj.vecPer2);

      elseif strcmpi(typeA, 'optical')

        obj.A = @(x) FEPack.tools.opticalPotential(x, obj.dualVec1, obj.dualVec2, true);

      elseif strcmpi(typeA, 'trigonometric')

        obj.A = @(x) obj.trigonometricPolynomial(x, true);

      elseif strcmpi(typeA, 'custom')

        if (nargin < 3)
          error(['Option', typeA, ' was selected: a function handle has to be provided.']);
        end

        obj.A = funA;

      elseif strcmpi(typeA, 'none')

        obj.A = [];
      
      else

        error(['Option ', typeA, ' unrecognized; only options are ''atomic'', ''optical'', ''trigonometric'', ''custom''']);

      end

      % W
      if strcmpi(typeW, 'optical')

        obj.W = @(x) FEPack.tools.opticalPotential(x, obj.dualVec1, obj.dualVec2, false);

      elseif strcmpi(typeW, 'trigonometric')

        obj.W = @(x) obj.trigonometricPolynomial(x, false);

      elseif strcmpi(typeW, 'custom')

        if (nargin < 4)
          error(['Option', typeW, ' was selected: a function handle has to be provided.']);
        end

        obj.W = funW;

      elseif strcmpi(typeW, 'none')

        obj.W = [];
      
      else

        error(['Option ', typeW, ' unrecognized; only options are ''atomic'', ''optical'', ''trigonometric'', ''custom''']);

      end

    end
    
    function val = trigonometricPolynomial(obj, x, is_even)
      
      if (is_even)
        val = cos(x(:, 1:2) *  obj.dualVec1) +...
              cos(x(:, 1:2) *  obj.dualVec2) +...
              cos(x(:, 1:2) * (obj.dualVec1  + obj.dualVec2));
      else
        val = sin(x(:, 1:2) *  obj.dualVec1) +...
              sin(x(:, 1:2) *  obj.dualVec2) +...
              sin(x(:, 1:2) * (obj.dualVec1  + obj.dualVec2));
      end

    end

  end

end
