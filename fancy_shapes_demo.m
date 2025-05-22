%% Compile of experiments with synthetic data. 
% Run the code and we will be asked to choose one of the synthetic datasets in the paper. 
% The ambient dimension is set to D = 100 (line 20), sigma is set to .05 (line 21) 
% and the overall noise, sqrt(3)*.05, is evenly distributed across the D dimensions. 
% The LAPD parameters are set in line 87-96. See main.m for descri
% If you see either the estimation of d or tau is very wrong, go to lines
% 149-158 to provide the corresponding oracle values. 

clear
close all
addpath(genpath(pwd));
rng("default")

prompt=['Please choose one fancy shape from below (in string format): \n ' ...
    'dollar sign, olympic rings, three curves, rose circles \n ' ...
    'two 2spheres, two triangles, three planes, swiss roll \n' ...
    'two 3spheres, two 4spheres, two 5spheres \n'];
DATAopts.shape = input(prompt); 
DATAopts.noise_type = 'uniform';   % Noise type. 
DATAopts.ambdim = 100;             % Ambient dimension. 
DATAopts.sigma = 0.05;             % Sigma; tau = sqrt(3)*sigma. 

% Set parameters for different shapes. 
if strcmp(DATAopts.shape, "dollar sign")
    DATAopts.number  = [5300, 2000]; DATAopts.intrdim = 1; 

elseif strcmp(DATAopts.shape, "olympic rings")
    DATAopts.number = [1200, 1200, 1200, 1200, 1200]; DATAopts.intrdim = 1;  

elseif strcmp(DATAopts.shape, "three curves")
    DATAopts.number = [5000, 5000, 5000]; DATAopts.intrdim = 1; 
    DATAopts.angles = [0, pi/3, 2*pi/3];   

elseif strcmp(DATAopts.shape, "rose circles")
    DATAopts.number = [3500, 1250, 1250]; DATAopts.intrdim = 1;  
    
elseif strcmp(DATAopts.shape, "two 2spheres")
    DATAopts.number = [3000, 3000];  DATAopts.intrdim = 2;  

elseif strcmp(DATAopts.shape, "two triangles")
    DATAopts.number = [6000, 6000];  DATAopts.intrdim = 2; 
    DATAopts.angles = [0, pi/10];   

elseif strcmp(DATAopts.shape, "three planes")
    DATAopts.number = [2000, 2000, 2000]; DATAopts.intrdim = 2;   
    DATAopts.angles = [0, pi/2, 4*pi/5];  

elseif strcmp(DATAopts.shape, "cone plane")
    DATAopts.number = [4000, 2000]; DATAopts.intrdim = 2;   

elseif strcmp(DATAopts.shape, "swiss roll")
    DATAopts.number = [4000, 2000]; DATAopts.intrdim = 2;    

elseif strcmp(DATAopts.shape, "two 3spheres")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 3;    

elseif strcmp(DATAopts.shape, "two 4spheres")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 4;    

elseif strcmp(DATAopts.shape, "two 5spheres")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 5;    

else
    error('Unknown shape.')
end

% Set LAPD parameters. 
LAPDopts.K = length(DATAopts.number); 
LAPDopts.intrdim = DATAopts.intrdim; 
LAPDopts.noise_level = sqrt(3)*DATAopts.sigma;
if strcmp(DATAopts.shape, "swiss roll") || strcmp(DATAopts.shape, "two triangles")
    LAPDopts.intrdim = DATAopts.intrdim; 
elseif strcmp(DATAopts.shape, "three curves")
    LAPDopts.noise_level = sqrt(3) * DATAopts.sigma; 
end
%LAPDopts.epsilon = 0.200;

% Generate random seeds. 
seeds = randsample(100, 20);

for i = 1:length(seeds)
    if i==1, disp('##########################################'); end
    if i==1, disp(strcat("Now testing ",DATAopts.shape," with tau = ", num2str(sqrt(3)*DATAopts.sigma), " using ", num2str(length(seeds)), " seeds.")); end
    DATAopts.rngSeed = seeds(i);
    
    % Generate data    
    [data, labelsGT] = simdata(DATAopts);

    % Run LAPD
    if exist("LAPDopts","var")
        [intrinsicDim,epsilon,k_hat,labels,time,misc] = main(data,LAPDopts);
    else
        [intrinsicDim,epsilon,k_hat,labels,time,misc] = main(data);
    end

    % Calculate accuracy
    [OA]= accuracy(labels, labelsGT);
    disp(strcat("Clustering accuracy for seed No.", num2str(i) ," = ", num2str(OA), "."))
    
    record(i,:) = [OA time epsilon misc.denoisingcutoff];
    if  ~exist("LAPDopts.intrdim", "var") || ~exist("LAPDopts.sigma","var")
        disp("_______________________________________________")
    end
    if i==length(seeds), disp('###############################################'); end
end

[mean(record(:,1)) std(record(:,1)) mean(record(:,2))]
if DATAopts.intrdim == 1 || DATAopts.intrdim == 2
    comparison(data,DATAopts,labels,labelsGT)
end