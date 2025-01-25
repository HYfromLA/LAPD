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
    DATAopts.number = [500, 500]; DATAopts.intrdim = 2;

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


% Set LAPD parameters. 
LAPDopts.K = length(DATAopts.number); 
LAPDopts.intrdim = DATAopts.intrdim; 
LAPDopts.noise_level = DATAopts.noise_level;
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

    % Run PBC
    [labels, ~,time] = Path_Based_Clustering(data,325,100,45,LAPDopts.K); 

    % Calculate accuracy
    [OA]= GetAccuracies(labels, labelsGT, LAPDopts.K);
    disp(strcat('Clustering accuracy for seed No.', num2str(i) ,' = ', num2str(OA), '.'))
    record(i,:) = [OA time];
    if  ~exist("LAPDopts.intrdim", "var") || ~exist("LAPDopts.noise_level","var")
        disp("_______________________________________________")
    end
    if i==20, disp('###############################################'); end
end

[mean(record(:,1)) std(record(:,1)) mean(record(:,2))]

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

%1d,pi/2:0:100,100,10 .01:200,100,50 .02:375,100,50 .03:550,100,50
%.04:750,100,50 .05:925,100,50 .06:1200,100,50 .07:1200,100,100
%1d,pi/4:0:100,100,10 .01:250,300,50 .02:500,300,50 .03:700,300,50
%.04:700,300,50 .05:700,300,50
%1d,pi/6:0:100,100,10 .01:200,100,50
%2d,pi/2:0:150,100,10 .01:400,100,10 .02:300,100,20 .03:500,100,25
%.04:600,100,30 .05:700,100,35 .06:800,100,35 .07:950,100,40
%.08:950,100,40 %.09:950,100,40
%2d,pi/4:0:100,100,10 .01:100,100,20 .02:200,100,25 .03:400,100,25
%.04:500,100,30 .05:600,100,35 .06:? 
%2d,pi/6:0:100,100,10 .01:200,100,10 .02:300,100,25 .04:410,100,30
%.04:400,100,35 .05:450,100,40 .06:500,100,45

%DS:225,100,45; OR:90,140,50; TC:500,150,50; RC:50,150,50; TS:85,300,30; TT:200,100,45; TP:275,100,45 SR: 3S:450,100,45; 4S:750,150,45; 5S:750,150,50;

load("mnist-full.mat")
classid = ismember(labelsGT,[0:6]); % Pick subset of all the classes. 
X=X(classid,:); labelsGT=labelsGT(classid,:);  
[labels, ~,time] = Path_Based_Clustering(X,4000,100,30,7); 
OA = accuracy(labels, labelsGT); 
[OA time]

% Coil20
%[0:6,8:17,19] 275,700,90 .38644
%[0:17,19] 275,700,90 .40509
%[0:19] 325,700,90 .43681

% mnist_test
%[0:5] 7500,100,30 .452
%[0:7] 9500,100,30 .414
%[0:8] 10500,100,30 .387
%[0:9] 11500,100,30 .358

% mnist_full
%[0:6] 1100,100,30 .2588
%[0:8] ?,?,? ?
%[0:9] ?,?,? ?

% USPS
%[0:5] 1300,100,80 .440
%[0:6,8:9] 100,100,80 .424
%[0:9] 700,100,80, .378