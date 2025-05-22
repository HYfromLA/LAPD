function [I_new, J_new, W_new] = filter_edge_list(I, J, W, rowsTokeep)
% Filter and reindex edge list using a logical vector rowsTokeep
% Inputs:
%   I, J, W - original edge list (symmetric assumed)
%   rowsTokeep - logical vector of length nn (true for nodes to keep)
% Outputs:
%   I_new, J_new, W_new - filtered edge list with reindexed node labels

% Step 1: Keep only edges between kept nodes
keep_mask = rowsTokeep(I) & rowsTokeep(J);
I_kept = I(keep_mask);
J_kept = J(keep_mask);
W_kept = W(keep_mask);

% Step 2: Remap old indices to compressed indices
% Get mapping: original indices → new indices
old_to_new = zeros(length(rowsTokeep), 1);
old_to_new(rowsTokeep) = 1:nnz(rowsTokeep);

% Step 3: Apply remapping
I_new = old_to_new(I_kept);
J_new = old_to_new(J_kept);
W_new = W_kept;
end