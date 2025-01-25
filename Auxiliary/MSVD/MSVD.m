function [d,tau] = MSVD(X,parallel)
    tX = X';     
    EstDimopts = struct('NumberOfTrials',3,'verbose',0,'MAXDIM',size(X,2),'MAXAMBDIM',size(X,2),'Ptwise',0,'PtIdxs',0,'NetsOpts',[],'UseSmoothedS',0,'EnlargeScales',0,'Deltas',[],'KeepV',0,'kNN',[],'AllPtsFast',0,'DownSample',1,'RandomizedSVDThres',inf);
    [d,noise_est] = EstDim_MSVD(tX, parallel, EstDimopts);
    tau = sqrt(size(X,2)-d)*noise_est; 
end