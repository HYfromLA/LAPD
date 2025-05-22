function [sortedCCmatrix, th] = connectedcomponents(I,J,W,simplices,numscales)
  
% This function constructs the connected components matrix (CCmatrix) using all the raw simplices.  
% Inputs:  I: row indices of nonzero entries in the adjacency matrix. 
%          J: column indices of nonzero entries in the adjacency matrix. 
%          W: weights of nonzero entries in the adjacency matrix. 
%          Simplices: the set of valid simplices. 
%          NumScales: number of scales for building the CCmatrix. 
% Outputs: sortedCCmatrix: the original CCmatrix (a matrix containing the connected components at each scale).
%          th: original scales coming with the sortedCCmatrix. 

N = size(simplices, 1);

%% Determine thresholds.
tmin = prctile(W, 1); %Take 1st percentile as smallest scale
tmax = 1.00001*max(W); 
th = exp(linspace(log(tmin),log(tmax),numscales)); 
th=[0,th];

%% Threshold graph.
[CCmatrix] = AlternateConnComp(N,I,J,W,th);

%% Create Sorted Matrix of Component Indices
CCmatrix = cat(2, CCmatrix, simplices); 
sortedCCmatrix = sortrows(CCmatrix, numscales:-1:1);

end


function [CCs] = AlternateConnComp(N,I,J,W,th)
   % Replacement for graphconncomp.m.  A is an n x n adjacency matrix, corresponding
   % to a nearest neighbor graog.  The function identifies the S
   % connected components C. This is done ala (Pothen, A. and Fan, C.J., 
   % 1990. Computing the block triangular form of a sparse matrix. ACM 
   % Transactions on Mathematical Software (TOMS), 16(4), pp.303-324).  This
   % is based on code from the gptoolbox by Alec Jacobson.
   %
   %
   % Input:
   %   -A:  n x n adjacency matrix
   % Outputs:
   %   -S:  scalar number of connected components
   %   -C:  matrix of connected component labels
 
   numscales = length(th)-1; 
   CCs=zeros(N,numscales); 
 
   for i=2:(numscales+1)
       Idx = W <= th(i); R=I(Idx); C=J(Idx); E=W(Idx);
       G = sparse(R, C, E, N, N); G = max(G, G');
       [p,~,r] = dmperm(G+speye(size(G)));
       C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
       C(p) = C; CCs(:,i-1) = C';
   end
 end