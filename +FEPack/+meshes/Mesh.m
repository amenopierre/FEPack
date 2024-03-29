%> @file Mesh.m
%> @brief Contains the meshes.Mesh class.
% =========================================================================== %
%> @brief Abstract class for all meshes
%>
%> Each mesh consists of nodes (points in dimension 1, 2, or 3) which
%> are connected by elements.
%>
%> The abstract Mesh class provides the interface for this
%> functionality by storing the dimension of the nodes, as well as their
%> coordinates.
%>
%> Specific informations regarding the elements can be found in the
%> subclasses.
% =========================================================================== %
classdef Mesh < FEPack.FEPackObject
  % FEPack.meshes.Mesh < FEPack.FEPackObject

  properties (SetAccess = protected)

    %> @brief Dimension of the geometry
    dimension = 0;

    %> @brief Tells if a mesh is structured or not
    is_structured = 0;

    %> @brief 3x1 vector that contains the number of nodes per side.
    %> The i-th element is the number of nodes along the i-th direction
    %> (set to 0 if the mesh dimension is lower)
    numEdgeNodes = 0;

    %> @brief Total number of vertices in the mesh
    numPoints = 0;

    %> @brief Matrix that contains the coordinates of the points. The k-th
    %> point is represented by the k-th line
    points = [];

    %> @brief Number of segments
    %> NOTE: a segment might be a volumic or surfacic entity depending on
    %> the dimension
    numSegments = 0;

    %> @brief Matrix that contains the definition of the segments. The k-th
    %> segment is represented by the k-th line. Each line contains the indices
    %> of the nodes that form the segment
    segments = [];

    %> @brief Number of triangles
    %> NOTE: a triangle might be a volumic or surfacic entity depending on
    %> the dimension
    numTriangles = 0;

    %> @brief Matrix that contains the definition of the triangle. The k-th
    %> triangle is represented by the k-th line. Each line contains the indices
    %> of its vertices
    triangles = [];

    %> @brief Number of tetrahedra
    numTetrahedra = 0;

    %> @brief Matrix that contains the definition of the tetrahedra. The k-th
    %> tetrahedra is represented by the k-th line. Each line contains the
    %> indices of its vertices
    tetrahedra = [];

    %> @brief Number of corners (2^dimension)
    numCorners = 0;

    %> @brief Indices of the corners
    corners = [];

    %> @brief List of domains
    domains = [];

    %> @brief References for the mesh points. Allows one to isolate the edge
    %> points, the interior points, and the corners
    refPoints = [];

    %> @brief References for edge segments
    refSegments = [];

    %> @brief References for face triangles
    refTriangles = [];

    %> @brief References for tetrahedra
    refTetrahedra = [];

    %> @brief As the meshes are periodic, these vectors map the boundary points
    %> along the x, y, z boundaries
    maps = {};

    %> @brief Relevant for cuboids, this vector maps the references of 
    %> the domains facing each other
    mapdomains =  [];

    % %> @brief Lagrange FE order
    % FEorder = [];

    % %> @brief Total number of degreees of freedom
    % numDOFglo = [];

  end

  methods (Abstract)

  end

  methods

    %> @brief domain() returns a FEDomain object given its name
    function dom = domain(mesh, domain_name)

      % List of names
      numDom = length(mesh.domains);
      names = cell(numDom, 1);
      for idom = 1:numDom
        names{idom} = mesh.domains{idom}.name;
      end

      % Find the id of the domain from its name
      domId = find(strcmp(names, domain_name));
      if (isempty(domId))
        error(['Le nom ''', domain_name, ''' ne figure pas dans la liste des domaines.']);
      end

      % Return the corresponding domain
      dom = mesh.domains{domId};

    end

    % %> @brief Choose the order of Lagrange FE
    % function attachFEorder(mesh, order)

    %   if (order < 1 || order > 2)
    %     error('L''ordre des éléments finis doit être compris entre 1 et 2.');
    %   end

    %   mesh.FEorder = order;

    %   switch (mesh.dimension)
    %   case 1
    %     mesh.numDOFglo = mesh.numPoints + (order - 1)*mesh.numSegments;
    %   case 2
    %     mesh.numDOFglo = mesh.numPoints +  (order - 1)*mesh.numSegments + ((order - 1)*(order - 2)/2)*mesh.numTriangles;
    %   case 3
    %     mesh.numDOFglo = mesh.numPoints +  (order - 1)*mesh.numSegments + ((order - 1)*(order - 2)/2)*mesh.numTriangles + ((order - 1)*(order - 2)*(order - 3)/6)*mesh.numTetrahedra;
    %   otherwise
    %     error('La dimension du maillage doit être comprise entre 1 et 3.');
    %   end

    % end

  end

end
