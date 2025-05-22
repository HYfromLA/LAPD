function [d,epsilon,k_hat,labels,time,misc] = main(X, LAPDopts)

% The main function of LAPD. 
% Inputs:  X: n * D n the number of samples; D the number of features. 
%          Options: d: the intrinsic dimension of the manifolds, if available. 
%                   K: the number of clusters, if known a priori. 
%                   noise: the noise level, i.e. if gaussian N(0, noise*I).
%                   epsilon: the inner radius of the annulus KNN graph. 
%                   bandwidth: the KNN bandwidth. 
%                   numScales: the number of scales for CCmatrix, default: 100. 
%                   denoisingmethod: 'ExamineFigure', 'Automatic_Elbow', or 'Automatic_Connectedness'. Default: 'Automatic_Elbow'.
%                   clusteringMethod: either 'Dendrogram' or 'Spectral', default: 'Dendrogram'.
%                   componentsize: a threshold to denoise very tiny component, default: 0.01. 
%                   allconnected: if 1, remove isolated components and keep only the giant graph. If 0 some components may not be connected. Default: 1. 
%                   parallel: 1 if parallel computing and 0 otherwise. 
% Outputs: d: the (estimated) intrinsic dimension of the manifolds. 
%          epsilon: the inner radius of the annulus KNN graph.
%          k_hat: (estimated) number of clusters. 
%          labels: cluster labels assigned to the data points. 
%          time: the runtime. 
%          misc: numsimplices: number of valid simplices. 
%                percentkept: percent of valid simplices surviving the quality filter. 
%                k1s: nearest neighbor used for each node. 
%                k2s: farthest neighbor used for each node. (k2s = k1s+bandwidth-1).
%                WLAPD: the largest within LAPD. 
%                BLAPD: the smallest between LAPD. 
%                SKNN: simplices' K-th nearest neighbor for denoising. 
%                denoisingcutoff: the denoising threshold. 


addpath(genpath(pwd)); rng("default")

if exist("LAPDopts",'var')
    if isfield(LAPDopts, 'intrdim'),          d = LAPDopts.intrdim; end 
    if isfield(LAPDopts, 'K'),                K = LAPDopts.K; end
    if isfield(LAPDopts, 'weight'),           weight = LAPDopts.weight; end
    if isfield(LAPDopts, 'noise_level'),      tau = LAPDopts.noise_level; end
    if isfield(LAPDopts, 'epsilon'),          epsilon = LAPDopts.epsilon; end 
    if isfield(LAPDopts, 'bandwidth'),        bandwidth = LAPDopts.bandwidth; end 
    if isfield(LAPDopts, 'filter'),           filter = LAPDopts.filter; end 
    if isfield(LAPDopts, 'numscales'),        numscales = LAPDopts.numscales; end
    if isfield(LAPDopts, 'denoisingmethod'),  denoisingmethod = LAPDopts.denoisingmethod; end
    if isfield(LAPDopts, 'SKNN'),             SKNN = LAPDopts.SKNN; end
    if isfield(LAPDopts, 'componentsize'),    componentsize = LAPDopts.componentsize; end
    if isfield(LAPDopts, 'parallel'),         parallel = LAPDopts.parallel; end
end

%% Setting default values to some parameters. 
if ~exist('weight', 'var'),            weight = "one sided"; end
if ~exist('denoisingmethod', 'var'),   denoisingmethod = 'automatic'; end %'automatic_elbow', 'examinefigure'
if ~exist('componentsize', 'var'),     componentsize = 0.01; end
if ~exist('parallel', 'var'),          parallel = 0; end
if parallel, parpool; end
tic; 

%% Estimating d, tau, and computing default epsilon. 

if ~exist('d','var') && ~exist('tau','var') 
    [d,tau] = MSVD(X,parallel); 
    fprintf('The intrinsic dimension found = %d \n', d);
    fprintf('The noise level found = %.5f. \n', tau);
elseif ~exist('d','var') && exist('tau','var') 
    [d,~] = MSVD(X,parallel); 
    fprintf('The intrinsic dimension found = %d \n', d);
elseif exist('d','var') && ~exist('tau','var') && ~exist('epsilon','var')
    [~,tau] = MSVD(X,parallel); 
    fprintf('The noise level found = %.5f. \n', tau);
end

[n, ~] = size(X); [IDXs,Dists] = knnsearch(X,X,'K',floor(0.1*n));

if ~exist('epsilon','var') 
    epsilon = sqrt(2)*tau;
    epsilon_L = prctile(Dists(:,ceil(2.5*log(n))), 90);
    samples = X(randsample(1:n,min(250,n)),:); epsilon_U = max(0.25*max(pdist(samples)),epsilon_L); 
    epsilon = min(max(epsilon, epsilon_L), epsilon_U); 
    clear samples 
end

if ~exist('bandwidth','var'), bandwidth = 25 ; end
if ~exist('numscales','var'), numscales = 50; end 

%% Building epsilon graph for each data point. 
[~, k1s]=max(Dists > epsilon, [], 2); k2s = k1s + bandwidth - 1;
while max(k2s) > n
    Prompt=['Warning: epsilon = ', num2str(epsilon),' is too large. Please enter a smaller epsilon... \n'];
    epsilon = input(Prompt); [~, k1s]=max(Dists > epsilon, [], 2); k2s = k1s + bandwidth - 1;
end
newIDXs = zeros(n, bandwidth); newDists = zeros(n, bandwidth); 
for j=1:n, newIDXs(j,:) = IDXs(j, k1s(j):k2s(j)); newDists(j,:) = Dists(j, k1s(j):k2s(j)); end
clear Dists k2s

%% Constructing simplices and calculate LAPD between adjacent simplices. 
if exist("filter","var" )
    [simplices,sharedfaces,posinsim,percentkept] = buildsimplex(X,n,d,newIDXs,newDists,parallel,filter);
else
    [simplices,sharedfaces,posinsim,percentkept] = buildsimplex(X,n,d,newIDXs,newDists,parallel);
end
nn = size(simplices,1);
clear newIDXs newDists
[I,J,W] = adjacency(X,d,weight,simplices,sharedfaces,posinsim,parallel);
clear X sharedfaces posinsim

%% Constructing CCmatrix. 
if ~exist("SKNN","var"), SKNN=floor(10*log(n)); end
[sortedCCmatrix, th] = connectedcomponents(I,J,W,simplices,numscales);
clear simplices
[denoisedsimplices, rowsTokeep, cutoff] = denoise(n,d,sortedCCmatrix,th,denoisingmethod,SKNN);
clear sortedCCmatrix th
[I, J, W] = filter_edge_list(I, J, W, rowsTokeep);
[denoisedCCmatrix, th] = connectedcomponents(I,J,W,denoisedsimplices,numscales);
clear I J W

%% Cluster
if exist("K","var")
    [k_hat,labels,latest_scale] = cluster(denoisedCCmatrix,numscales,n,componentsize,IDXs,K);
else
    [k_hat,labels,latest_scale] = cluster(denoisedCCmatrix,numscales,n,componentsize,IDXs);
end
clear IDXs

misc.numsimplices = nn; misc.epsilon = epsilon; misc.k1s = k1s; misc.th = th;
misc.percentkept = percentkept; misc.SKNN = SKNN; misc.denoisingcutoff = cutoff;

time = toc; 
if parallel, delete(gcp); end

end