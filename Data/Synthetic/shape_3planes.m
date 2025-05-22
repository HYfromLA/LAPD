function[data,labels] = shape_3planes(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
rotation = DATAopts.angles;
tau = DATAopts.sigma; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed); 

labels = zeros(sum(n),1);
% For each plane
for j = 1:length(n)

  % Generate rotation
  c = cos(rotation(j));
  s = sin(rotation(j));
  G = eye(D);
  G(1:(d+1),1:(d+1)) = [ c  0  s;
                         0  1  0;
                        -s  0  c];

  % Generate noiseless, standardized data
  thisData = zeros([n(j) D]);
  %thisData(:,1:d) = -0.5 + rand([n(j) 2]);
   
  if j==3
      temp = -0.5 + rand([n(j) 2]); temp(:,1) = 1.5*temp(:,1);
      thisData(:,1:d) = temp; %1.5*(-0.5 + rand([n(j) 2]));
  elseif j == 1
      temp = -0.5 + rand([n(j) 2]); temp(:,1) = 1.25*temp(:,1);
      thisData(:,1:d) = temp; %-0.5 + rand([n(j) 2]);
  else
      temp = -0.5 + rand([n(j) 2]); temp(:,1) = 1.25*temp(:,1);
      thisData(:,1:d) = temp; %-0.5 + rand([n(j) 2]);
  end

  % Pollute standardized data with noise
  switch noise_type
    case 'gaussian'
      thisData = thisData + normrnd(0,tau/sqrt(D-d),size(thisData));    %noise_stdev/sqrt(D-d)*randn(size(thisData));
    case 'uniform'
      thisData = thisData + sqrt(3)*tau/sqrt(D-d)*(2*rand(size(thisData))-1);
    otherwise
      error("Unknown noise type")
  end

  % Rotate polluted data
  thisData = thisData * G;

  % Re-locate plane 2 and 3. 
  if j==2 
      thisData(:,1) = thisData(:,1) -0.300;  %-0.175;
      thisData(:,3) = thisData(:,3) +0.250;  %+0.205;
  end

  if j==3 
      thisData(:,1) = thisData(:,1) + 0.100; %+0.055
      thisData(:,3) = thisData(:,3) + 0.180; %+0.150
  end

  % Place data in right locations in output.
  row1 = 0 + sum(n(1:(j-1))) + 1;
  row2 = sum(n(1:j));
  data(row1:row2,:) = thisData;
  labels(row1:row2) = j;

  plot3(data(:,1),data(:,2),data(:,3),'.')

end