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


%options = struct('SKNN',floor(2*d*log(n))); % Default values for kappa

% Parse name-value pairs
%for i = 1:2:length(varargin)
%    if isfield(options, varargin{i})
%        options.(varargin{i}) = varargin{i+1};
%    elseif strcmp(varargin{i}, 'K')
%        K = varargin{i+1}; % Capture 'e' if provided
%    else
%        error('Invalid option: %s', varargin{i});
%    end
%end
%SKNN = options.SKNN;

if ~exist("SKNN","var"), SKNN=floor(5*d*log(n)); end

%% Construct CCmatrix. 

[sortedCCmatrix, th] = connectedcomponents(I,J,W,simplices,numscales);
clear simplices
[denoisedCCmatrix,newTh,cutoff] = denoise(n,d,I,J,W,sortedCCmatrix,th,denoisingmethod,SKNN);
clear sortedCCmatrix I J W

%{
if exist("K","var")
    [goodclustering,matching,K,labels_S,finalCCmatrix] = clusterscheck(d,denoisedCCmatrix,componentsize,K);  
    
    if goodclustering && ~matching         
        [goodclustering,matching,K_est,labels_S,finalCCmatrix] = clusterscheck(d,denoisedCCmatrix,componentsize);          
        disp(strcat("User requested ",num2str(K), " clusters but we can only find ",num2str(K_est), " clusters."))
        K = K_est; 
    end
else
    [goodclustering,matching,K,labels_S,finalCCmatrix] = clusterscheck(d,denoisedCCmatrix,componentsize);  
end
%}
end