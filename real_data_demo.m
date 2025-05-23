addpath(genpath(pwd));

%% COIL20. 
clear
close all
load("coil20.mat")
classid = ismember(labelsGT,[1:20]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);  
LAPDopts.intrdim = 1; LAPDopts.K = length(unique(labelsGT));
LAPDopts.epsilon = 0; %LAPDopts.denoisingmethod = 'examinefigure'; 
LAPDopts.bandwidth = 22; LAPDopts.weight = "two sided";
[intrinsicDim,epsilon,k_hat,labels,time,misc] = main(X,LAPDopts);
OA = accuracy(labels,labelsGT); 

%[1:7,9:18,20] bandwidth = 17, denoising=0.8: 98.15%
%[1:18,20] bandwidth = 19, denoising=1: 93.06%
%[1:20] bandwidth = 22, denoising=1: 90.42%

%[1:7,9:18,20] bandwidth = 20, denoising=default: 96.37%
%[1:18,20] bandwidth = 20, denoising=default: 91.30%
%[1:20] bandwidth = 22, denoising=default: 90.28%

%% USPS.
clear
close all
load("USPS.mat")
classid = ismember(labelsGT,[0:9]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);  
LAPDopts.intrdim = 2; LAPDopts.K = length(unique(labelsGT));
LAPDopts.epsilon = 0; LAPDopts.weight = "two sided"; LAPDopts.bandwidth = 17;
LAPDopts.filter = 1.3; LAPDopts.componentsize = 0.008; 
[intrinsicDim,epsilon,k_hat,labels,time,misc] = main(X,LAPDopts);
OA = accuracy(labels, labelsGT); 

% [0:6] bandwidth=17,denoising=default, 97.1%.    
% [0:6,8:9] bandwidth=17,LAPDopts.filter = 1.3, denoising=default, 94.2%.
% [0:9] bandwidth=17,LAPDopts.filter = 1.3,denoising=default,LAPDopts.componentsize = 0.008, 93.5%.

%% MNIST-test. 
clear
close all
load("mnist_test.mat")
classid = ismember(labelsGT,[0:9]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);
LAPDopts.intrdim = 2; LAPDopts.K = length(unique(labelsGT));
LAPDopts.epsilon = 0; %LAPDopts.denoisingmethod = "examinefigure";
LAPDopts.bandwidth = 14; LAPDopts.weight = "two sided";  
[intrinsicDim,epsilon,k_hat,labels,time,misc] = main(X,LAPDopts);
OA= accuracy(labels, labelsGT);

%[0:7] bandwidth = 14, denoising=default, SKNN=default: 96.13%
%[0:8] bandwidth = 14, denoising=default, SKNN=default: 94.32%
%[0:9] bandwidth = 14, denoising=default, SKNN=default: 80.98%

%% MNIST-full. 
clear
close all
load("mnist_full.mat")
classid = ismember(labelsGT,[0:6]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);
LAPDopts.intrdim = 2; LAPDopts.K = length(unique(labelsGT));
%LAPDopts.filter = 1.20; LAPDopts.numscales = 100; 
LAPDopts.epsilon = 0; %LAPDopts.denoisingmethod = 'examinefigure'; 
LAPDopts.bandwidth = 11; LAPDopts.weight = "two sided"; LAPDopts.knnnumber = 30;
[intrinsicDim,epsilon,k_hat,labels,time,misc] = main(X,LAPDopts);
OA= accuracy(labels, labelsGT);

%[0:6] bandwidth = 10, denoising=default: 93.87%   551.83
%[0:6,8] bandwidth = 10, 95.59%  669.88
%[0:8] bandwidth = 15,  filter=1.2, denoising=1.2, numscales=100 : .950
%[0:9] bandwidth = 15, filter=1.2, denoising=1.2, SKNN=default: .859