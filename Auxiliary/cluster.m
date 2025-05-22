function [n_component, node_labels, latest_scale] = cluster(CCmatrix,m,n,min_percentage,IDXs,n_component)
    % CLUSTER - Perform node clustering based on connected component matrix.
    % Inputs:
    %   CCmatrix - Connected component matrix of size (n, m + d + 1).
    %   m - Number of scales in CCmatrix.
    %   n - Total number of data points (nodes).
    %   n_component - Desired number of clusters (optional, provide as empty []).
    %   min_percentage - Minimum percentage threshold for valid labels.
    %
    % Outputs:
    %   node_labels - A vector of length n with cluster labels (0, 1, ...).
    %   latest_scale - The index of the column (scale) used for clustering.

    % Extract scales and simplices
    scales = CCmatrix(:, 1:m);

    label_counts_per_scale = zeros(m, 1);
    label_to_columns = containers.Map('KeyType', 'double', 'ValueType', 'any');

    % Count valid labels for each scale
    for scale = 1:m
        simplex_labels = scales(:, scale);
        [unique_labels, ~, label_idx] = unique(simplex_labels);
        counts = accumarray(label_idx, 1);
        percentages = counts / length(simplex_labels);

        % Filter valid labels
        valid_labels = unique_labels(percentages >= min_percentage);
        num_valid_labels = length(valid_labels);
        label_counts_per_scale(scale) = num_valid_labels;

        if ~isKey(label_to_columns, num_valid_labels)
            label_to_columns(num_valid_labels) = [];
        end
        %label_to_columns(num_valid_labels) = [label_to_columns(num_valid_labels), scale];
        current_columns = label_to_columns(num_valid_labels);
        if ~ismember(scale, current_columns) % Ensure no duplicate entries
            label_to_columns(num_valid_labels) = [current_columns, scale];
        end
    end

    if ~exist("n_component", "var")
    % Get all unique cluster sizes and their counts
    cluster_sizes = cell2mat(keys(label_to_columns));
    cluster_counts = cellfun(@length, values(label_to_columns));
   
    % Find the cluster size with the maximum count
    valid_indices = cluster_sizes > 1; 
    valid_sizes = cluster_sizes(valid_indices); 
    valid_sizes
    valid_counts = cluster_counts(valid_indices); 
    valid_counts

    if isempty(valid_sizes)
        error('No valid cluster size greater than 1.')
    end

    max_count = max(valid_counts);
    max_count_candidates = valid_sizes(valid_counts == max_count);
   
    % Choose the smallest cluster size among those with max count
    n_component = min(max_count_candidates);
   
    % If the chosen n_component is 1, find the smallest size > 1
    if n_component == 1
        valid_candidates = cluster_sizes(cluster_sizes > 1);
        if ~isempty(valid_candidates)
            n_component = min(valid_candidates);
        else
            error('No valid cluster size greater than 1 found.');
        end
    end
   
    fprintf('Estimated n_component: %d\n', n_component);
    end


    % Match or adjust n_component
    if isKey(label_to_columns, n_component)
        matched_columns = label_to_columns(n_component);
    else
        candidates = cell2mat(keys(label_to_columns));
        larger_candidates = candidates(candidates > n_component);
        smaller_candidates = candidates(candidates < n_component & candidates > 1);

        if ~isempty(larger_candidates) && ~isempty(smaller_candidates)
            closest_larger = min(larger_candidates, [], 'omitnan');
            larger_count = length(label_to_columns(closest_larger));

            closest_smaller = max(smaller_candidates, [], 'omitnan');
            smaller_count = length(label_to_columns(closest_smaller));

            if larger_count > smaller_count
                n_component = closest_larger;
            else
                n_component = closest_smaller;
            end
        elseif ~isempty(larger_candidates) && isempty(smaller_candidates)
            closest_larger = min(larger_candidates, [], 'omitnan');
            n_component = closest_larger; 
        elseif isempty(larger_candidates) && ~isempty(smaller_candidates)
            closest_smaller = max(smaller_candidates, [], 'omitnan');
            n_component = closest_smaller; 
        end

        %if larger_count > smaller_count
        %    n_component = closest_larger;
        %elseif smaller_count > larger_count || isempty(closest_larger)
        %    n_component = closest_smaller;
        %else
        %    n_component = closest_smaller;
        %end
        fprintf('Adjusted n_component to %d.\n', n_component);
        matched_columns = label_to_columns(n_component);
    end

    % Use the latest matching scale
    latest_scale = matched_columns(end);
    simplex_labels = scales(:, latest_scale);

    % Filter valid labels based on min_percentage
    valid_labels = unique(simplex_labels);
    percentages = histcounts(simplex_labels, [valid_labels; max(valid_labels)+1]) / length(simplex_labels);
    valid_labels = valid_labels(percentages >= min_percentage);

    % Filter CCmatrix for valid simplices
    is_valid = ismember(simplex_labels, valid_labels);
    CCmatrix_filtered = CCmatrix(is_valid, :);
    simplices_filtered = CCmatrix_filtered(:, m+1:end);
    simplex_labels_filtered = CCmatrix_filtered(:, latest_scale);

    % Relabel valid labels consecutively starting from 0
    [~, ~, new_labels] = unique(simplex_labels_filtered);
    %relabeled_simplices = new_labels - 1;

    % Initialize node labels
    node_labels = zeros(n, 1);

    % Assign labels to nodes using majority voting
    for node = 1:n
        simplex_indices = any(simplices_filtered == node, 2);
        simplex_labels_for_node = new_labels(simplex_indices);

        if ~isempty(simplex_labels_for_node)
            unique_labels = unique(simplex_labels_for_node);
            counts = histcounts(simplex_labels_for_node, [unique_labels; max(unique_labels)+1]);
            [~, max_idx] = max(counts);
            node_labels(node) = unique_labels(max_idx);
        end
    end
    
    % if any of the data are removed from the simplices set and have not
    % been labeled from majority voting. 
    if any(node_labels == 0) 
        removednodes = find(node_labels == 0); 

        for i = 1:length(removednodes)
            NNs = IDXs(removednodes(i), :); labeledneighbors = []; j = 0; 
            while length(labeledneighbors) < 11 && j<500
                j=j+1; 
                if node_labels(NNs(j)) > 0
                    labeledneighbors = [labeledneighbors; NNs(j)];
                end
            end
            node_labels(removednodes(i)) = mode(node_labels(labeledneighbors));
        end
    end
end
