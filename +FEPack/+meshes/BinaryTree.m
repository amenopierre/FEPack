%> @file BinaryTree.m
%> @brief Contains the meshes.BinaryTree class.
% =========================================================================== %
%> @brief Binary tree structure used to locally refine 1D intervals
%> based on splitting rules that may involve an indicator function
%>
%> Coded with the help of DeepSeek
% =========================================================================== %
classdef BinaryTree < FEPack.FEPackObject
  % FEPack.meshes.BinaryTree < FEPack.meshes.Mesh

  properties

    %> @brief Binary tree nodes [xmin, xmax, depth]
    leaves

    %> @brief vector containing values within each leaf
    leaf_vals

    %> @brief Maximum refinement depth
    maxDepth 

    %> @brief Splitting rules
    splitting_rules

    %> @brief Nx1 vector of x evaluation points
    cache_coords

    %> @brief Nx1 vector of function values
    cache_vals 

  end
  
  methods

    function obj = BinaryTree(maxDepth, rootInterval, rules)
      % function obj = BinaryTree(maxDepth, rootInterval)
      % Initialize binary tree with root node
      %
      % inputs: * maxDepth (integer), maximum depth
      %         * rootInterval 2-sized vector containing 
      %                [xmin, xmax, depth]
      %         * rules (function handle) that takes 
      %                 as inputs 
      %               - a Nxp vector representing
      %               values (p = 2 is the number of endpoints)
      %               - a Nx1 vector containing leaves depths
      %               and returns a Nx1 boolean 
      if (nargin < 3)
        % Default rule: force early split, then split leaves where
        % the indicator function is above a certain threshold
        threshold = 1e2;
        force_depth = min(2, maxDepth);
        blow_up_rule = @(fvals) (any(abs(fvals) > threshold, 2) &...
                                ~all(abs(fvals) > threshold, 2));
        early_split_rule = @(depth) (depth <= force_depth);
        rules = @(fvals, depth) early_split_rule(depth) | blow_up_rule(fvals);
      end
      if (nargin < 2)
        rootInterval = [0, 1];
      end

      obj.leaves = [rootInterval, 0];
      obj.maxDepth = maxDepth;
      obj.cache_coords = zeros(0, 1);
      obj.cache_vals = zeros(0, 1);
      obj.set_splitting_rules(rules);
    end
    
    function obj = record_value(obj, P, val)
      % function obj = RECORD_VALUE(obj, P, val)
      % Records function evaluation in vector storage
      %
      % INPUTS: * obj, FEPack.meshes.BinaryTree object.
      %         * P, a N-by-1 matrix containing the 
      %              coordinates of the points.
      %         * val, a N-by-1 vector containing 
      %                the corresponding values
      
      obj.cache_coords = [obj.cache_coords; P];
      obj.cache_vals   = [obj.cache_vals; val];
    end
    
    function [val, exists] = get_cached(obj, P, almostzero)
      % function [val, exists] = GET_CACHED(obj, P, almostzero)
      %
      % Check if coordinate is included in cache
      %
      % INPUTS: * obj, FEPack.meshes.BinaryTree object.
      %         * P, a N-by-1 matrix containing the 
      %              coordinates of the points.
      %         * almostzero (optional) tolerance.
      %
      % OUPUTS: * val, a Ncached-by-1 vector containing
      %           the cached values. Ncached <= N is the
      %           number of points that are already in
      %           the cached coordinates.
      %
      %         * exists, a N-by-1 vector. exists(i) is one
      %           if P(i) is in the cache, and 0 otherwise.         
      
      if (nargin < 3)
        almostzero = 1e-8;
      end

      Pc = obj.cache_coords;
      N  = size(P,  1);
      Nc = size(Pc, 1);

      distance = abs(P * ones(1, Nc) - ones(N, 1) * Pc');
      
      % Find points contained in the cache
      matches = (distance < almostzero);
      exists = any(matches, 2);
      
      % Modify value accordingly
      idExists = matches(exists, :)';
      Nex = size(idExists, 2);
      val = obj.cache_vals * ones(1, Nex);
      val = val(idExists);
    end
    
    function obj = refine(obj, indicator_function)
      % Main refinement loop with cached evaluations

      for iter = 1:obj.maxDepth
        newLeaves = [];
        nodes = obj.leaves;
        
        % Find leaves that should be split
        should_split = obj.check_node(nodes, indicator_function);

        % Nodes that are not split are added to the new leaves
        newLeaves = [newLeaves; nodes(~should_split, :)]; %#ok

        % Split leaves that should be
        newLeaves = [newLeaves; obj.split_node(nodes(should_split, :))]; %#ok

        obj.leaves = newLeaves;
        % obj.visualize();
        % pause;
      end

      obj.compute_leaf_values;
    end
    
    function set_splitting_rules(obj, rules)
      % obj = SET_SPLITTING_RULES(obj, rules)
      % Define the rules used for splitting
      %
      % inputs: * rules (function handle) that takes 
      %                 as inputs 
      %               - a Nxp vector representing
      %               values (p = 2 is the number of endpoints)
      %               - a Nx1 vector containing leaves depths
      %               and returns a Nx1 boolean 
      obj.splitting_rules = rules;
    end

    function should_split = check_node(obj, nodes, indicator_function)
      % Node evaluation with caching
      xmin  = nodes(:, 1);
      xmax  = nodes(:, 2);
      depth = nodes(:, 3);
      
      % Check endpoints
      endpoints{1} = xmin;
      endpoints{2} = xmax;

      N  = size(nodes, 1);
      fvals = zeros(N, 2);
      
      for idJ = 1:2
        % Check if function has been already evaluated at endpoints 
        [cached_val, exists] = obj.get_cached(endpoints{idJ});

        fvals( exists, idJ) = cached_val;
        fvals(~exists, idJ) = indicator_function(endpoints{idJ}(~exists, :));
        obj = obj.record_value(endpoints{idJ}(~exists, :), fvals(~exists, idJ));
      end
      
      should_split = obj.splitting_rules(fvals, depth);
    end
    
    function children = split_node(~, nodes)
      % Splits node into 2 children
      xmidpoints = (nodes(:, 1) + nodes(:, 2))/2;
      newDepth = nodes(:, 3) + 1;
      
      N = size(nodes, 1);
      children = zeros(2*N, 3);

      children(1:2:2*N, :) = [nodes(:,1), xmidpoints, newDepth];
      children(2:2:2*N, :) = [xmidpoints, nodes(:,2), newDepth];
    end
    
    function compute_leaf_values(obj)
      numLeaves = size(obj.leaves, 1);

      % Match cache points to leaves and compute mean values
      obj.leaf_vals = NaN(numLeaves, 1);

      for idI = 1:numLeaves
        node = obj.leaves(idI, :);
        
        % Find all cache points within this leaf
        in_leaf = obj.cache_coords(:, 1) >= node(1) & ...
                  obj.cache_coords(:, 1) <= node(2);
        
        if any(in_leaf)
          obj.leaf_vals(idI) = mean(obj.cache_vals(in_leaf));
        end
      end

      if all(isnan(obj.leaf_vals))
        error('No cached values found in binary tree leaves');
      end
    end

    function visualize(obj)
      axis equal;
      hold on;

      for idI = 1:size(obj.leaves, 1)
        plot([obj.leaves(idI, 1), obj.leaves(idI, 2)], [0, 0], 'ro');
      end
    end

    function visualize_cache(obj, fun)
      % Visualize cached values on the binary tree mesh as a 2D plot
      % Uses same rendering approach as visualize_function but with cached values
      if (nargin < 2)
          fun = @(x) x;
      end

      % Get valid leaves (non-NaN values)
      val_is_NaN = isnan(obj.leaf_vals);
      valid_leaves = obj.leaves(~val_is_NaN, :);
      valid_vals = fun(obj.leaf_vals(~val_is_NaN));
      
      % Get NaN leaves and assign value -1
      nan_leaves = obj.leaves(val_is_NaN, :);
      nan_vals = -1 * ones(size(nan_leaves, 1), 1);
      
      % Sort leaves by x-coordinate for proper plotting
      [~, sortIds] = sort(valid_leaves(:, 1));
      sorted_leaves = valid_leaves(sortIds, :);
      sorted_vals = valid_vals(sortIds);
      
      [~, sortNanIds] = sort(nan_leaves(:, 1));
      sorted_nan_leaves = nan_leaves(sortNanIds, :);
      sorted_nan_vals = nan_vals(sortNanIds);
      
      % Create step function representation for valid values
      X = zeros(2 * size(sorted_leaves, 1), 1);
      Y = zeros(2 * size(sorted_leaves, 1), 1);
      
      for idI = 1:size(sorted_leaves, 1)
          idJ = 2*(idI-1) + 1;
          X(idJ)   = sorted_leaves(idI, 1);
          X(idJ+1) = sorted_leaves(idI, 2);
          Y(idJ)   = sorted_vals(idI);
          Y(idJ+1) = sorted_vals(idI);
      end
      
      % Create step function representation for NaN values
      Xnan = zeros(2*size(sorted_nan_leaves, 1), 1);
      Ynan = zeros(2*size(sorted_nan_leaves, 1), 1);
      
      for idI = 1:size(sorted_nan_leaves, 1)
          idJ = 2*(idI - 1) + 1;
          Xnan(idJ)   = sorted_nan_leaves(idI, 1);
          Xnan(idJ+1) = sorted_nan_leaves(idI, 2);
          Ynan(idJ)   = sorted_nan_vals(idI);
          Ynan(idJ+1) = sorted_nan_vals(idI);
      end
      
      % Plot the step functions
      plot(X, Y, 'b-');
      hold on;
      plot(Xnan, Ynan, 'k-');
      
      xlabel('x');
      ylabel('Function Value');
      title('Cached Function Values on Binary Tree Mesh');
      grid on;
      
      % Add markers at the center of each interval for clarity
      centers = mean(sorted_leaves(:, 1:2), 2);
      nan_centers = mean(sorted_nan_leaves(:, 1:2), 2);
      
      plot(centers, sorted_vals, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
      plot(nan_centers, sorted_nan_vals, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
      
      % % Add legend
      % legend('Valid values', 'NaN values (set to -1)', 'Location', 'best');
      
      hold off;
    end

    function [XnanBB, Xpeaks] = find_nan_boundary_and_peaks(obj, threshold, cluster_threshold_ratio)

      if (nargin < 3)
        cluster_threshold_ratio = 1.5;
      end

      % Sort leaves by spatial coordinate
      [~, sortKey] = sort(obj.leaves(:, 1));
      sorted_leaves = obj.leaves(sortKey, :);
      sorted_vals = obj.leaf_vals(sortKey);

      % Find indices of NaN leaves
      nanBool = isnan(sorted_vals);

      % Start and end positions of consecutive blocks of NaN values
      XnanBB = [sorted_leaves(find(diff([0; nanBool]) ==  1), 1),...
                sorted_leaves(find(diff([nanBool; 0]) == -1), 2)];

      % Peaks indices
      vals = obj.leaf_vals;
      vals(isnan(vals)) = -Inf;

      peaksIds = find(vals >= threshold);
      XpeaksInit = (obj.leaves(peaksIds, 1) + obj.leaves(peaksIds, 2)) / 2;
      XpeaksVals = vals(peaksIds);

      % Group peaks indices into clusters, and take average in each cluster
      N = numel(XpeaksInit);

      if (N == 0)
        Xpeaks = [];
        return;
      elseif (N == 1)
        Xpeaks = XpeaksInit;
        return;
      end

      % Xdiff = diff(XpeaksInit);
      cluster_threshold = cluster_threshold_ratio * min(obj.leaves(:, 2) - obj.leaves(:, 1));
      
      clusterIds = ones(N, 1);
      currCluster = 1;

      for idI = 2:N
        if ((XpeaksInit(idI) - XpeaksInit(idI-1)) > cluster_threshold)
          currCluster = currCluster + 1;
        end

        clusterIds(idI) = currCluster;
      end

      numClusters = clusterIds(end);
      Xpeaks = zeros(numClusters, 1);

      for idI = 1:numClusters
        mask = (clusterIds == idI);
        Xclusters = XpeaksInit(mask);
        % [~, maxId] = max(XpeaksVals(mask));
        % Xpeaks(idI) = Xclusters(maxId);
        Xpeaks(idI) = mean(Xclusters);
      end
    end

  end
end