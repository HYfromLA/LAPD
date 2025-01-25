function [I,J,W] = adjacency(X,d,weight,simplices,sharedfaces,posinsim,parallel)

% This function construct the adjacency matrix using the angle distance min(theta, pi-theta). 
% Inputs:  X: the n * D datastes. 
%          d: the intrinsic dimension of the manifolds. 
%          simplices: the set of valid simplices. 
%          sharedfaces: all the shared faces by the simplices. 
%          posinsim: the position of all the simplices in Simplices sharing a particular shared face. 
%          parallel: 1 if parallel computing and 0 otherwise. 
% Outputs: I: row indices of nonzero entries in the adjacency matrix. 
%          J: column indices of nonzero entries in the adjacency matrix. 
%          W: weights of nonzero entries in the adjacency matrix. 
% Note: the adjacency can be easily built as A = sparse(I, J, W, N, N); A = max(A, A');
% We don't build it until when it is necessary in order to save memory. 

nn = size(sharedfaces, 1); 

if parallel
    parfor i = 1:nn
        sharednodes = sharedfaces(i, :); vecrows = posinsim{i}; m=length(vecrows);
        block = simplices(vecrows, :); block = block';
        othernodes = setdiff(block, sharednodes,'stable'); uniquenodes = cat(1, sharednodes', othernodes); 
        coords = X(uniquenodes,:); temp = (coords - coords(1,:))'; temp = temp(:, 2:end);

        [V1, V2] = pairvectors(d, m, temp);

        tempmat = repmat((1:m)', 1, m); lower = tril(tempmat, -1); 
        lower = lower(:); lower = lower(lower>0); cols = vecrows(lower); 
        cumsum = linspace(0, (m-1)*m, m); repcumsum = repelem(cumsum, (m-1):-1:0)'; 
        saverep2 = lower+repcumsum; V2 = V2(saverep2,:);
        thetas = real(acos(dot(V1,V2,2))); Thetas{i} = thetas; 
        rows = repelem(VecRows(1:(m-1)), (m-1):-1:1, 1);
        I{i} = rows; J{i} = cols; 
    end
else
    for i = 1:nn
        sharednodes = sharedfaces(i, :); vecrows = posinsim{i}; m=length(vecrows); 
        block = simplices(vecrows, :); block = block';
        othernodes = setdiff(block, sharednodes,'stable'); uniquenodes = cat(1, sharednodes', othernodes); 
        coords = X(uniquenodes,:); temp = (coords - coords(1,:))'; temp = temp(:, 2:end);

        [V1, V2] = pairvectors(d, m, temp);

        tempmat = repmat((1:m)', 1, m); lower = tril(tempmat, -1); 
        lower = lower(:); lower = lower(lower>0); 
        cumsum = linspace(0, (m-1)*m, m); repcumsum = repelem(cumsum, (m-1):-1:0)'; 
        saverep2 = lower+repcumsum; V2 = V2(saverep2,:);
        thetas = real(acos(dot(V1,V2,2))); Thetas{i} = thetas; 
        rows = repelem(vecrows(1:(m-1)), (m-1):-1:1, 1); cols = vecrows(lower); 

        I{i} = rows; J{i} = cols; 
    end
end

I = cat(1, I{:}); J = cat(1, J{:}); Thetas = cat(1, Thetas{:}); 
if strcmp(weight, "two sided")
    W = min(pi-Thetas,Thetas); W(W < 1e-8) = 1e-8;
else
    W = pi-Thetas; W(W < 1e-8) = 1e-8; idx = W<(pi/2)*1.001; W=W(idx); I=I(idx); J=J(idx); 
end

end


function [V1, V2] = pairvectors(d, m, temp)

% This function pairs adjacent simplices and normalizes their feature vectors. 
% Inputs:  d: intrinsic dimension of manifolds. 
%          m: the number of simplices that share the shared face. 
%          temp: 
% Outputs: V1,V2: normalized feature vectors of all simplices that share the same shared face. 

if d == 1
    NormedVecs = normc(temp)';
    V1 = repelem(NormedVecs(1:(m-1), :), (m-1):-1:1, 1);
    V2 = repmat(NormedVecs, (m-1), 1);
else
    SharedVectors = temp(:,1:(d-1)); OtherVectors = temp(:, d:end);   % The shared vectors and the other vectors. 
    % Make shared vectors an orthonormal basis; use residuals from
    % projecting other vectors to the basis to calculate the angles. 
    P = orth(SharedVectors); Q = P';  
    Residuals = OtherVectors - P*(Q*OtherVectors); Residuals = Residuals ./ vecnorm(Residuals, 2, 1);    %Residuals = normc(OtherVectors - P*(Q*OtherVectors));
    Residuals = Residuals';
    V1 = repelem(Residuals(1:(m-1), :), (m-1):-1:1, 1);
    V2 = repmat(Residuals, (m-1), 1);
end

end