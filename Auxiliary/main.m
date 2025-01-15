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


tic; addpath(genpath(pwd)); rng("default")
X = double(X); 

if exist("LAPDopts",'var')
    if isfield(LAPDopts, 'intrdim'),          d = LAPDopts.intrdim; end 
    if isfield(LAPDopts, 'K'),                K = LAPDopts.K; end
    if isfield(LAPDopts, 'weight'),           weight = LAPDopts.weight; end
    if isfield(LAPDopts, 'noise_level'),      noise = LAPDopts.noise_level; end
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
if ~exist('denoisingmethod', 'var'),   denoisingmethod = 'automatic_elbow'; end %'automatic_connectedness', 'automatic_elbow', 'examinefigure'
if ~exist('componentsize', 'var'),     componentsize = 0.01; end

%% Estimating d, tau, and computing default epsilon. 
if ~exist('d','var') 
    tX = X';     
    EstDimopts = struct('NumberOfTrials',3,'verbose',0,'MAXDIM',size(X,2),'MAXAMBDIM',size(X,2),'Ptwise',0,'PtIdxs',0,'NetsOpts',[],'UseSmoothedS',0,'EnlargeScales',0,'Deltas',[],'KeepV',0,'kNN',[],'AllPtsFast',0,'DownSample',1,'RandomizedSVDThres',inf);
    fprintf('Intrinsic dimension is not provided. Now estimating... \n');
    [d,noise_est] = EstDim_MSVD(tX, parallel, EstDimopts);
    fprintf('The intrinsic dimension found is %d \n', d);
    if ~exist('noise','var') 
        noise = sqrt(size(X,2)-d)*noise_est; 
        fprintf('The noise level found is %.5f. \n', noise);
    end
    clear tX
end

if ~exist('epsilon','var') % if only epsilon is not given by user    
    if ~exist('noise','var') % estimate noise
        tX = X';     
        EstDimopts = struct('NumberOfTrials',3,'verbose',0,'MAXDIM',size(X,2),'MAXAMBDIM',size(X,2),'Ptwise',0,'PtIdxs',0,'NetsOpts',[],'UseSmoothedS',0,'EnlargeScales',0,'Deltas',[],'KeepV',0,'kNN',[],'AllPtsFast',0,'DownSample',1,'RandomizedSVDThres',inf);
        fprintf('Noise level is not provided. Now estimating... \n');
        [~,noise_est] = EstDim_MSVD(tX, parallel, EstDimopts);
        noise = sqrt(size(X,2)-d)*noise_est; 
        fprintf('The noise level found is %.5f. \n', noise);
        clear tX
    end
    
    epsilon = 2.5*d^(0.5)*noise;  %sqrt(d)
    [n, ~] = size(X); [IDXs,Dists] = knnsearch(X,X,'K',floor(0.5*n));
    epsilon_L = mean(Dists(:,max(21, floor(log(n))))); 
    samples = X(randsample(1:n,min(250,n)),:); epsilon_U = max(0.25*max(pdist(samples)),epsilon_L); 
    epsilon = min(max(epsilon, epsilon_L), epsilon_U); 
    clear samples 
end

if ~exist('bandwidth','var'), bandwidth = 20+5*(d-1); end

epsilon
if ~exist('numscales','var'),numscales = 100;end 

%% Building epsilon graph for each data point. 
if ~exist('Dists','var'), [n, ~] = size(X); [IDXs,Dists] = knnsearch(X,X,'K',floor(0.5*n)); end
[~, k1s]=max(Dists > epsilon, [], 2); k2s = k1s + bandwidth - 1;
mean(k1s)
while max(k2s) > n
    Prompt=['Warning: epsilon = ', num2str(epsilon),' is too large. Please enter a smaller epsilon... \n'];
    epsilon = input(Prompt);
    [~, k1s]=max(Dists > epsilon, [], 2); k2s = k1s + bandwidth - 1;
end
newIDXs = zeros(n, bandwidth); newDists = zeros(n, bandwidth); 
for j=1:n, newIDXs(j,:) = IDXs(j, k1s(j):k2s(j)); newDists(j,:) = Dists(j, k1s(j):k2s(j)); end
clear Dists k2s

%% Constructing simplices and calculate LAPD between adjacent simplices. 
if exist("filter","var" )
    [simplices,sharedfaces,posinsim,percentkept] = buildsimplex(X,n,d,newIDXs,newDists,filter);
else
    [simplices,sharedfaces,posinsim,percentkept] = buildsimplex(X,n,d,newIDXs,newDists);
end
nn = size(simplices,1);
clear newIDXs newDists
[I,J,W] = adjacency(X,d,weight,simplices,sharedfaces,posinsim);
clear X sharedfaces posinsim

%% Constructing CCmatrix. 
if exist("SKNN","var")
    [denoisedCCmatrix,newTh,SKNN,cutoff] = buildCCmatrix(d,n,I,J,W,simplices,numscales,denoisingmethod,SKNN); 
else
    [denoisedCCmatrix,newTh,SKNN,cutoff] = buildCCmatrix(d,n,I,J,W,simplices,numscales,denoisingmethod); 
end
clear I J W simplices 

%% Cluster
if exist("K","var")
    [k_hat,labels,latest_scale] = cluster(denoisedCCmatrix,numscales,n,componentsize,IDXs,K);
else
    [k_hat,labels,latest_scale] = cluster(denoisedCCmatrix,numscales,n,componentsize,IDXs);
end
clear IDXs





%if exist("K","var")
%    if exist("SKNN","var")
%        [goodclustering,k_hat,denoisedCCmatrix,newTh,cutoff,SKNN,labels_S] = buildCCmatrix(d,n,I,J,W,simplices,numscales,denoisingmethod,componentsize,parallel,"SKNN",SKNN,"K",K);
%    else
%        [goodclustering,k_hat,denoisedCCmatrix,newTh,cutoff,SKNN,labels_S] = buildCCmatrix(d,n,I,J,W,simplices,numscales,denoisingmethod,componentsize,parallel,"K",K);
%    end
%else
%    if exist("SKNN","var")
%        [goodclustering,k_hat,denoisedCCmatrix,newTh,cutoff,SKNN,labels_S] = buildCCmatrix(d,n,I,J,W,simplices,numscales,denoisingmethod,componentsize,parallel,"SKNN",SKNN);
%    else
%        [goodclustering,k_hat,denoisedCCmatrix,newTh,cutoff,SKNN,labels_S] = buildCCmatrix(d,n,I,J,W,simplices,numscales,denoisingmethod,componentsize,parallel);
%    end
%end
%clear I J W simplices

%% Performing clustering. 
%if goodclustering
%   [labels,WLAPDs,WLAPD,BLAPDs,BLAPD] = label_dendrogram(d,n,denoisedCCmatrix,newTh,labels_S,IDXs);
%else
%    disp('Current parameter setting cannot partition the data.')
%    k_hat = 1; labels = ones(n,1); WLAPD = NaN; BLAPD = NaN; 
%end

misc.numsimplices = nn; 
misc.epsilon = epsilon; 
misc.percentkept = percentkept;
misc.k1s = k1s; 
%misc.WLAPDs = WLAPDs; misc.WLAPD = WLAPD; 
%misc.BLAPDs = BLAPDs; misc.BLAPD = BLAPD; 
misc.SKNN = SKNN; 
misc.denoisingcutoff = cutoff; %misc.dnsCCmatrix = denoisedCCmatrix; 

time = toc; 
%if parallel, delete(gcp); end

end