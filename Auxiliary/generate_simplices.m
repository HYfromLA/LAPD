function [simplices, shared_edges, edge_to_simplices] = generate_simplices(X, e, q, d)
% X: n x D data matrix
% e: lower bound for edge length
% q: contraction factor (< 1), upper bound is e/q
% d: intrinsic dimension
% Output:
% simplices: m x (d+1) matrix, each row is a simplex
% shared_edges: list of shared edges (pairs of node indices)
% edge_to_simplices: cell array, each entry is the indices of simplices sharing that edge

n = size(X, 1);
% Step 1: Compute pairwise distances and find valid edges using pdist2
Dists = triu(pdist2(X, X),1); % n x n distance matrix
[i_idx, j_idx] = find(Dists >= e & Dists <= e/q); % upper triangle only to avoid duplicates
valid_edges = [i_idx, j_idx];
clear Dists i_idx j_idx

% Map valid edges for fast lookup
edge_map = containers.Map;
for k = 1:size(valid_edges,1)
    key = edge_key(valid_edges(k,1), valid_edges(k,2));
    edge_map(key) = true;
end

% Step 2: Generate all (d+1)-simplices from combinations of nodes
simplices_list = [];
combs = nchoosek(1:n, d+1);

for i = 1:size(combs,1)
    nodes = combs(i,:);
    all_valid = true;
    pairs = nchoosek(nodes, 2);
    for j = 1:size(pairs,1)
        key = edge_key(pairs(j,1), pairs(j,2));
        if ~isKey(edge_map, key)
            all_valid = false;
            break;
        end
    end
    if all_valid
        simplices_list = [simplices_list; sort(nodes)];
    end
end

simplices = simplices_list;
m = size(simplices,1);

% Step 3: Identify shared edges and corresponding simplex indices
edge_index_map = containers.Map;
for i = 1:m
    simplex = simplices(i,:);
    pairs = nchoosek(simplex, 2);
    for j = 1:size(pairs,1)
        key = edge_key(pairs(j,1), pairs(j,2));
        if isKey(edge_index_map, key)
            edge_index_map(key) = [edge_index_map(key), i];
        else
            edge_index_map(key) = i;
        end
    end
end

% Collect shared edges and their associated simplices
shared_edges = [];
edge_to_simplices = {};
keys_list = keys(edge_index_map);
for i = 1:length(keys_list)
    key = keys_list{i};
    simplex_ids = edge_index_map(key);
    if numel(simplex_ids) > 1
        ids = sscanf(key, '%d_%d');
        shared_edges = [shared_edges; ids'];
        edge_to_simplices{end+1} = simplex_ids;
    end
end

end

function k = edge_key(a, b)
% Create a unique string key for an undirected edge
a_b = sort([a, b]);
k = sprintf('%d_%d', a_b(1), a_b(2));
end
