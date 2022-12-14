%> @file CellBVP.m
%> @brief Contains the +solver.CellBVP class.
% =========================================================================== %
%> @brief class for solver of boundary value problem in a cell
% =========================================================================== %
classdef CellBVP < FEPack.solver.BVPObject
  % FEPack.solver.CellBVP < FEPack.solver.BVPObject

  properties (SetAccess = protected)

    %> @brief Mesh of the domain
    mesh = [];

    %> @brief Bilinear integral
    bilinearIntg = [];

    %> @brief Linear integral
    linearIntg = [];

    %> @brief Essential conditions
    ecs = [];

    %> @brief Solution
    vec = [];

  end

  methods

    % Initialization of the solver
    function initialize(cellBVPsolver, dimension, varargin)

      % Preliminary check on the dimension
      validDimension = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 4);
      if ~validDimension(dimension)
          error('La dimension doit être un scalaire entier valant 1, 2, ou 3.');
      end

      % Default values
      defaultNumNodes = (2^(7-dimension)) * ones(1, dimension);
      u = FEPack.pdes.PDEObject;
      v = dual(u);
      defaultBilinearIntegrand = (grad(u)*grad(v)) + id(u)*id(v);
      defaultLinearIntegrand = id(v);

      % Argument validation
      validNumNodes = @(x) isnumeric(x) && isvector(x) && (length(x) == dimension);

      % Add inputs
      entrees = inputParser;
      entrees.addRequired('dimension', validDimension)
      entrees.addParameter('numEdgeNodes', defaultNumNodes, validNumNodes);
      entrees.addParameter('BilinearIntegrand', defaultBilinearIntegrand, @(x) isa(x, 'FEPack.pdes.Form'));
      entrees.addParameter('LinearIntegrand', defaultLinearIntegrand, @(x) isa(x, 'FEPack.pdes.Form'));
      

      % Parse the inputs
      parse(entrees, dimension, varargin{:});
      cellBVPsolver.mesh = entrees;

    end

    % Resolution
    function solve(cellBVPsolver)






    end
    % Representation of the output

  end

end
