%% COIL20. 
clear
close all
load("coil20.mat")
classid = ismember(labelsGT,[1:20]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);  
LAPDopts.intrdim = 1; LAPDopts.K = length(unique(labelsGT));
LAPDopts.epsilon = 0; LAPDopts.denoisingmethod = 'examinefigure'; 
LAPDopts.bandwidth = 22; LAPDopts.weight = "two sided"; 
[intrinsicDim,epsilon,k_hat,labels,time,misc] = main(X,LAPDopts);
OA = accuracy(labels,labelsGT); 

%[1:7,9:18,20] bandwidth = 17, denoising=0.8: 98.23%
%[1:18,20] bandwidth = 19, denoising=1: 92.84%
%[1:20] bandwidth = 22, denoising=1: 90.42%

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
% [0:6,8:9] bandwidth=17,denoising=default,LAPDopts.filter = 1.3, 94.2%.
% [0:9] bandwidth=17,denoising=default,LAPDopts.filter = 1.3,LAPDopts.componentsize = 0.008, 93.5%.

%% MNIST-test. 
clear
close all
load("mnist_test.mat")
classid = ismember(labelsGT,[0:9]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);
LAPDopts.intrdim = 2; LAPDopts.K = length(unique(labelsGT));
LAPDopts.epsilon = 0; LAPDopts.denoisingmethod = 'examinefigure';
%LAPDopts.filter = 1.2;
LAPDopts.bandwidth = 15; LAPDopts.weight = "two sided";  
[intrinsicDim,epsilon,k_hat,labels,time,misc] = main(X,LAPDopts);
OA= accuracy(labels, labelsGT);

%[0:7] bandwidth = 15, denoising=1.15, SKNN=default: 95.80%
%[0:8] bandwidth = 15, denoising=1.15, SKNN=default: 93.71%  %.9226
%[0:9] bandwidth = 15, denoising=1.15, SKNN=default: 80.98%

%% MNIST-full. 
clear
close all
load("mnist_test.mat")
classid = ismember(labelsGT,[0:9]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);
LAPDopts.intrdim = 2; LAPDopts.K = length(unique(labelsGT));
LAPDopts.filter = 1.2; LAPDopts.epsilon = 0;
%LAPDopts.denoisingmethod = 'examinefigure'; 
LAPDopts.bandwidth = 15; LAPDopts.weight = "two sided";  
[intrinsicDim,epsilon,k_hat,labels,time,misc] = main(X,LAPDopts);
OA= accuracy(labels, labelsGT);

%[0:5] bandwidth = 18, denoising=1.1, SKNN=default: .978
%[0:6] bandwidth = 18, denoising=1.1, SKNN=default: .965
%[0:7] bandwidth = 15, filter=1.2, denoising=1.2, SKNN=default: .971
%[0:8] bandwidth = 15, filter=1.2, denoising=1.2, SKNN=default: .957
%[0:9] bandwidth = 15, filter=1.2, denoising=1.2, SKNN=default: .859

% One sided
%[0:5] bandwidth = 20, filter=1.2, denoising=1.5, SKNN=default: .950
%[0:6] bandwidth = 20, filter=1.2, denoising=1.5, SKNN=default: .948
%[0:7] bandwidth = 18, filter=1.2, denoising=1.5, SKNN=default: .615