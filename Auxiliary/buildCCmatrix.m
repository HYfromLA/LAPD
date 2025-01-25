function [denoisedCCmatrix,newTh,SKNN,cutoff] = buildCCmatrix(d,n,I,J,W,simplices,numscales,denoisingmethod,SKNN)

% This function tries different number of scales to make sure that the built CCmatrix is clusterable. 
% If K is given, we also make sure that denoisedCCNums contains K if provided. 
% Inputs:  d: intrinsic dimension of manifolds. 
%          n: dataset size. 
%          I: row indices of nonzero entries in the adjacency matrix. 
%          J: column indices of nonzero entries in the adjacency matrix. 
%          W: weights of nonzero entries in the adjacency matrix. 
%          simplices: all the valid simplices.  
%          numscales: number of scales for building the CCmatrix. 
%          denoisingmethod: examinefigure, automatic_elbow, or automatic_connectedness. 
%          componentsize: component size threshold below which the components are deemed tiny and noisy. 
%          parallel: Upper end for adjusting the scale. 
%          K: user-requested number of clusters. 
% Outputs: goodclustering: the denoised CCmatrix. 
%          matching: new scales coming with the Denoisedmatrix. 
%          denoisedCCmatrix: unique knn distances to k2-th nearest neighbors. 
%          newTh: cutoff to remove simplices with large knn distances. 
%          k2.  

if ~exist("SKNN","var"), SKNN=floor(5*d*log(n)); end
[sortedCCmatrix, th] = connectedcomponents(I,J,W,simplices,numscales);
clear simplices
[denoisedCCmatrix,newTh,cutoff] = denoise(n,d,I,J,W,sortedCCmatrix,th,denoisingmethod,SKNN);
clear sortedCCmatrix I J W
end