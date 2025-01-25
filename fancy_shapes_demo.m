%% Compile of experiments with synthetic data.  
% If you see either the estimation of d or tau is very wrong, go to lines
% 149-158 to provide the corresponding oracle values. 

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
DATAopts.ambdim = 100;

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
    DATAopts.number = [6000, 6000]; %5000,5000
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
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 2;

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

elseif strcmp(DATAopts.shape, "two cuboids")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 3;

    prompt='Please choose an intersection angle in [pi/6, pi/2] to test: \n ';
    angle = input(prompt); 
    DATAopts.angles = [0, angle];

    if angle == pi/2
        max_noise = 0.08; 
    elseif angle == pi/3
        max_noise = 0.08; 
    elseif angle == pi/4
        max_noise = 0.05; 
    elseif angle == pi/5
        max_noise = 0.04; 
    elseif angle == pi/6
        max_noise = 0.02; 
    end

elseif strcmp(DATAopts.shape, "two tesseracts")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 4;

    prompt='Please choose an intersection angle in [pi/6, pi/2] to test: \n ';
    angle = input(prompt); 
    DATAopts.angles = [0, angle];

    if angle == pi/2
        max_noise = 0.07; 
    elseif angle == pi/3
        max_noise = 0.08; 
    elseif angle == pi/4
        max_noise = 0.05; 
    elseif angle == pi/5
        max_noise = 0.04; 
    elseif angle == pi/6
        max_noise = 0.02; 
    end

elseif strcmp(DATAopts.shape, "two penteracts")
    DATAopts.number = [3000, 3000]; DATAopts.intrdim = 5;

    prompt='Please choose an intersection angle in [pi/6, pi/2] to test: \n ';
    angle = input(prompt); 
    DATAopts.angles = [0, angle];

    if angle == pi/2
        max_noise = 0.06; 
    elseif angle == pi/3
        max_noise = 0.08; 
    elseif angle == pi/4
        max_noise = 0.05; 
    elseif angle == pi/5
        max_noise = 0.04; 
    elseif angle == pi/6
        max_noise = 0.02; 
    end

else
    error('Unknown shape.')
end

prompt=['Please choose a noise level in [0, ', num2str(max_noise),']: \n '];
DATAopts.noise_level = input(prompt); 


% Set LAPD parameters. 
LAPDopts.K = length(DATAopts.number); 
LAPDopts.intrdim = DATAopts.intrdim; 
LAPDopts.noise_level = DATAopts.noise_level;
LAPDopts.denoisingmethod = 'automatic_elbow'; 
%LAPDopts.epsilon = 0.23;
%LAPDopts.bandwidth = 20;
LAPDopts.weight = "two sided";
if strcmp(DATAopts.shape, "swiss roll") || strcmp(DATAopts.shape, "two triangles")
    LAPDopts.intrdim = DATAopts.intrdim; 
elseif strcmp(DATAopts.shape, "three curves")
    LAPDopts.noise_level = DATAopts.noise_level; 
end

% Generate random seeds. 
seeds = randsample(100, 20);

for i = 1:length(seeds)
    if i==1, disp('##########################################'); end
    if i==1, disp(strcat("Now testing ",DATAopts.shape," with tau = ", num2str(DATAopts.noise_level), " using ", num2str(length(seeds)), " seeds.")); end
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

    % Confusion matrix. 
    %cm = confusionchart(labelsGT, labels);
    %cm.Title = 'Confusion Matrix';
    %cm.RowSummary = 'row-normalized'; % Show row-normalized values
    %cm.ColumnSummary = 'column-normalized'; % Show column-normalized values
    
    record(i,:) = [OA time epsilon misc.denoisingcutoff];
    if  ~exist("LAPDopts.intrdim", "var") || ~exist("LAPDopts.noise_level","var")
        disp("_______________________________________________")
    end
    if i==length(seeds), disp('###############################################'); end
end

[mean(record(:,1)) std(record(:,1)) mean(record(:,2))]

% Save the result as .mat files. 
%if strcmp(DATAopts.noise_type, "gaussian")
%    save(strcat("result_",DATAopts.shape,"_G.mat"),"record")
%elseif strcmp(DATAopts.noise_type, "uniform")
%    save(strcat("result_",DATAopts.shape,"_U.mat"),"record")
%end 


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