function [simplices,sharedfaces,posinsim,percentkept] = buildsimplex(X, n, d, IDXs, Dists, filter)

% This function construct valid simplices. 
% Inputs:  X: the n * D datastes. 
%          d: the intrinsic dimension of the manifolds. 
%          IDXs: the n * B matrix containing the epsilon nearest neighbors with bandwidth B. 
%          Dists: corresponding knn distances to IDXs. 
%          Parallel: 1 if parallel computing and 0 otherwise. 
% Outputs: Simplices: the set of valid simplices. 
%          SharedFaces: all the shared faces by the simplices. 
%          PosInSim: the position of all the simplices in Simplices sharing a particular shared face. 
%          PercentKept: the percent of simplices surviving the filter. 100% for d=1.  

parallel = 0;

bandwidth = size(IDXs, 2); %filter = 1+0.1*d; %2:1.20;3:1.25;4:   %1.15+0.15*(d-1);


if ~exist("filter", "var")
    switch d
        case 2, filter = 1.25; 
        case 3, filter = 1.33;
        case 4, filter = 1.41;
        case 5, filter = 1.49;
    end
end

if d == 1
    IDXs = IDXs'; simplices = IDXs(:); sharedfaces = (1:n)'; 
    simplices = cat(2, repelem(sharedfaces, bandwidth, 1), simplices);
    [simplices, ~, ~] = unique(sort(simplices,2),'rows');  
    percentkept = 1; %largest_edgelength = vecnorm(X(simplices(:,2),:)-X(simplices(:,1),:),2,2);
else
    nodesperm = nchoosek(1:d, 2); numperm = nchoosek(d, 2); % Each edge vector has two nodes. %NodesPerm2 = nchoosek(1:(d+1), d);
    m=nchoosek(bandwidth, d); 
    if parallel 
        parfor i = 1:n
            NNs = IDXs(i,:); DDs = Dists(i,:);
            nodes = nchoosek(NNs, d); node1edgesnorms = nchoosek(DDs, d);
            simplicestemp = cat(2, repelem(i, m)', nodes);
            node1s = nodes(:, nodesperm(:,1)); node1s = node1s(:);
            node2s = nodes(:, nodesperm(:,2)); node2s = node2s(:);
            thirdedgesnorms = vecnorm(X(node2s,:)-X(node1s,:),2,2);
            thirdedgesnorms = reshape(thirdedgesnorms, [], numperm); edgesnorms = cat(2, node1edgesnorms, thirdedgesnorms); % Each row contains the norms of all the edges.
            quality = max(edgesnorms, [], 2) ./ min(edgesnorms, [], 2); 
            goodquality = quality <= filter; simplices{i} = simplicestemp(goodquality, :); 
        end
    else
        for i = 1:n
            NNs = IDXs(i,:); DDs = Dists(i,:);
            nodes = nchoosek(NNs, d); node1edgesnorms = nchoosek(DDs, d);
            simplicestemp = cat(2, repelem(i, m)', nodes);
            node1s = nodes(:, nodesperm(:,1)); node1s = node1s(:);
            node2s = nodes(:, nodesperm(:,2)); node2s = node2s(:);

            %num_nodes = numel(node1s); num_chunks = ceil(num_nodes / chunk_size);          
            thirdedgesnorms = vecnorm(X(node2s,:)-X(node1s,:),2,2);
            thirdedgesnorms = reshape(thirdedgesnorms, [], numperm); edgesnorms = cat(2, node1edgesnorms, thirdedgesnorms); % Each row contains the norms of all the edges.
            max_edgelength = max(edgesnorms, [], 2); min_edgelength = min(edgesnorms, [], 2);
            goodquality = (max_edgelength ./ min_edgelength) <= filter; simplices{i} = simplicestemp(goodquality, :); 
            clear nodes node1edgesnorms simplicestemp node1s node2s thirdedgesnorms edgesnorms
        end
    end
    simplices = cat(1, simplices{:});

    % Keep only unique simplices. 
    simplices = unique(sort(simplices, 2), 'rows'); 
    percentkept = size(simplices, 1)/(n*nchoosek(bandwidth, d));
end

nodesperm2 = nchoosek(1:(d+1), d); shared = []; 
v = simplices(:, nodesperm2');  positions = repmat((1:size(simplices, 1))', (d+1), 1);
for j = 1 : d 
    w = v(:, j:d:end);  w = w(:); shared = cat(2,shared, w);
end
clear w v

[~, ~, T] = unique(shared, 'rows');
[a, ia] = sort(T); nn = length(a); positions = positions(ia);
clear T
groupEnd = [find(diff(a)==1); nn];  groupStart = [1; groupEnd(1:end-1)+1];
clear a
groups = cat(2, groupStart, groupEnd); 
save = groups(:, 1) ~= groups(:, 2);  % Keep only shared faces that appear more than once (so that at least an angle can form).

MoreThanOneAppearance = groups(save, :); 
clear groups
shared = shared(ia,:); sharedfaces = shared(MoreThanOneAppearance(:,1), :);

for i = 1:size(MoreThanOneAppearance, 1)
    posinsim{i} = sort(positions(MoreThanOneAppearance(i, 1):MoreThanOneAppearance(i, 2)));
end

end