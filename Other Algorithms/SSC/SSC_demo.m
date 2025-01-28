clear
close all
addpath(genpath(pwd));
rng("default")

prompt=['Please choose one fancy shape from below (in string format): \n ' ...
    'two lines, two planes, two cuboids, two tesseracts, two penteracts \n' ...
    'dollar sign, olympic rings, three curves, rose circles \n ' ...
    'two 2spheres, two triangles, three planes, swiss roll \n' ...
    'two 3spheres, two 4spheres, two 5spheres \n'];
DATAopts.shape = input(prompt); 

% Set the noise type. 
prompt=['Please choose one noise type from below (in string format): \n ' ...
    'uniform, gaussian \n'];
DATAopts.noise_type = input(prompt); 

% Set the ambient dimension. 
DATAopts.ambdim = 2500;

% Set parameters for different shapes. 
if strcmp(DATAopts.shape, "dollar sign")
    DATAopts.number  = [5300, 2000]; 
    DATAopts.intrdim = 1; 
    max_noise = 0.150; 

elseif strcmp(DATAopts.shape, "olympic rings")
    DATAopts.number = [1200, 1200, 1200, 1200, 1200]; 
    DATAopts.intrdim = 1;  
    max_noise = 0.070; 

elseif strcmp(DATAopts.shape, "three curves")
    DATAopts.number = [5000, 5000, 5000]; 
    DATAopts.angles = [0, pi/3, 2*pi/3];  
    DATAopts.intrdim = 1;  
    max_noise = 0.150; 

elseif strcmp(DATAopts.shape, "rose circles")
    DATAopts.number = [4000, 1000, 1000]; 
    DATAopts.intrdim = 1;  
    max_noise = 0.050; 
    
elseif strcmp(DATAopts.shape, "two 2spheres")
    DATAopts.number = [3000, 3000];  
    DATAopts.intrdim = 2;  
    max_noise = 0.140; 

elseif strcmp(DATAopts.shape, "two triangles")
    DATAopts.number = [5000, 5000]; 
    DATAopts.angles = [0, pi/10];   
    DATAopts.intrdim = 2; 
    max_noise = 0.05; 

elseif strcmp(DATAopts.shape, "three planes")
    DATAopts.number = [2000, 2000, 2000]; 
    DATAopts.angles = [0, pi/2, 4*pi/5];  
    DATAopts.intrdim = 2; 
    max_noise = 0.040; 

elseif strcmp(DATAopts.shape, "swiss roll")
    DATAopts.number = [4000, 2000]; DATAopts.intrdim = 2;    
    max_noise = 0.070; 

elseif strcmp(DATAopts.shape, "two 3spheres")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 3;    
    max_noise = 0.150; 

elseif strcmp(DATAopts.shape, "two 4spheres")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 4;    
    max_noise = 0.150; 

elseif strcmp(DATAopts.shape, "two 5spheres")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 5;    
    max_noise = 0.150; 

elseif strcmp(DATAopts.shape, "two lines")
    DATAopts.number = [2500, 2500]; DATAopts.intrdim = 1;

    prompt='Please choose an intersection angle in [pi/6, pi/2] to test: \n ';
    angle = input(prompt); 
    DATAopts.angles = [0, angle];

    if angle == pi/2
        max_noise = 0.13; 
    elseif angle == pi/3
        max_noise = 0.08; 
    elseif angle == pi/4
        max_noise = 0.05; 
    elseif angle == pi/5
        max_noise = 0.04; 
    elseif angle == pi/6
        max_noise = 0.04; 
    end

elseif strcmp(DATAopts.shape, "two planes")
    DATAopts.number = [2500, 2500]; DATAopts.intrdim = 2;

    prompt='Please choose an intersection angle in [pi/6, pi/2] to test: \n ';
    angle = input(prompt); 
    DATAopts.angles = [0, angle];

    if angle == pi/2
        max_noise = 0.15; %??
    elseif angle == pi/3
        max_noise = 0.08; %??
    elseif angle == pi/4
        max_noise = 0.05; %??
    elseif angle == pi/5
        max_noise = 0.04; %??
    elseif angle == pi/6
        max_noise = 0.02; %??
    end

    LAPDopts.intrdim = 2; 

elseif strcmp(DATAopts.shape, "two cuboids")
    DATAopts.number = [2500, 2500]; DATAopts.intrdim = 3;

    prompt='Please choose an intersection angle in [pi/6, pi/2] to test: \n ';
    angle = input(prompt); 
    DATAopts.angles = [0, angle];

    if angle == pi/2
        max_noise = 0.13; 
    elseif angle == pi/3
        max_noise = 0.08; 
    elseif angle == pi/4
        max_noise = 0.05; 
    elseif angle == pi/5
        max_noise = 0.04; 
    elseif angle == pi/6
        max_noise = 0.02; 
    end

    LAPDopts.intrdim = 3; 

else
    error('Unknown shape.')
end

prompt=['Please choose a noise level in [0, ', num2str(max_noise),']: \n '];
DATAopts.noise_level = input(prompt); 

% Generate random seeds. 
seeds = randsample(100, 20);

for i = 1:length(seeds)
    if i==1, disp('##########################################'); end
    if i==1, disp(strcat("Now testing ",DATAopts.shape," with tau = ", num2str(DATAopts.noise_level), " using ", num2str(length(seeds)), " seeds.")); end
    DATAopts.rngSeed = seeds(i);
    
    % Generate data    
    [data, labelsGT] = simdata(DATAopts);

    % Run SSC
    [missrate,labels, ~,time] = SSC(data',3,false,8000,false,1,labelsGT); 

    % Calculate accuracy
    [OA]= 1-missrate; 
    disp(strcat('Clustering accuracy for seed No.', num2str(i) ,' = ', num2str(OA), '.'))
    record(i,:) = [OA time];
    if  ~exist("LAPDopts.intrdim", "var") || ~exist("LAPDopts.noise_level","var")
        disp("_______________________________________________")
    end
    if i==length(seeds), disp('###############################################'); end
end

mean(record(:,2))

% Save the result as .mat files. 
if strcmp(DATAopts.noise_type, "gaussian")
    save(strcat("result_",DATAopts.shape,"_G.mat"),"record")
elseif strcmp(DATAopts.noise_type, "uniform")
    save(strcat("result_",DATAopts.shape,"_U.mat"),"record")
end 


%% Plots for comparison

viz.markersize = 15;

figure; 
label_ids = unique(labelsGT);
if (DATAopts.intrdim == 1 && strcmp(DATAopts.shape, "three curves")) || DATAopts.intrdim == 2
    subplot(1,2,1)
    hold on;
    for label = label_ids(:).'
        mask = (labelsGT == label);
        plot3(data(mask,1), data(mask,2), data(mask,3), '.', viz);
    end
    title('Oracle labels')
    hold off

    subplot(1,2,2)
    hold on;
    for label = label_ids(:).'
        mask = (labels == label);
        plot3(data(mask,1), data(mask,2), data(mask,3),'.', viz);
    end
    title('LAPD labels')
    hold off

elseif DATAopts.intrdim == 1
    subplot(1,2,1)
    hold on;
    for label = label_ids(:).'
        mask = (labelsGT == label);
        plot(data(mask,1), data(mask,2), '.', viz);
    end
    title('Oracle labels')
    hold off

    subplot(1,2,2)
    hold on;
    for label = label_ids(:).'
        mask = (labels == label);
        plot(data(mask,1), data(mask,2), '.', viz);
    end
    title('LAPD labels')
    hold off
end

%2lines 3,8000,1; 2planes 3,8000,1; 2cuboids 10,8000,1


%% COIL20
load("coil20.mat")
classid = ismember(labelsGT,[1:7,9:18,20]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);  
[~, ~, new_labels] = unique(labelsGT);
[missrate,labels, ~,time] = SSC(X',0,false,11,false,0.9,new_labels); 
1-missrate
% [1:7,9:18,20]: 0,9,0.7 .876; 0,11,0.7 .876;
% [1:18,20]: 0,7,0.7 .863; 

load("USPS.mat")
classid = ismember(labelsGT,[0:9]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);  
[~, ~, new_labels] = unique(labelsGT);
[missrate,labels, ~,time] = SSC(X',0,false,5,false,1,new_labels); 
accuracy(labels, labelsGT)

%[0:5] 5,1 [0:6,8:9] 5,1 [0:9] 

load("mnist_test.mat")
classid = ismember(labelsGT,[0:8]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);  
[~, ~, new_labels] = unique(labelsGT);
[missrate,labels, ~,time] = SSC(X',0,false,4,false,1,new_labels); 
accuracy(labels, labelsGT)

%[0:5] 4,1 [0:7] 4,1 %[0:8] 5,1 %[0:9] 5,1