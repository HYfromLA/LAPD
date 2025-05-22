function[data,labels] = shape_3curves(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
rotation = DATAopts.angles;
tau = DATAopts.sigma;  
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed);

%labels = zeros([sum(n) 1]);

% For each curve
for j = 1:length(n)

  % Generate rotation
  c = cos(rotation(j));
  s = sin(rotation(j));
  G = eye(D);
  G(1:2,1:2) = [ c  s; 
                -s  c];

  % Generate noiseless, standardized data
  thisData = randn(n(j), d+2); thisData(:,2) = 0;
  thisData = thisData ./ vecnorm(thisData,2,2);
  thisData(thisData(:,3)>0,:)=[]; nn(j)=size(thisData,1);
  thisData(:,d+3:D) = zeros(nn(j),D-d-2);

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

  % Place data in right locations in output.
  row1 = 0 + sum(nn(1:(j-1))) + 1;
  row2 = sum(nn(1:j));
  data(row1:row2,:) = thisData;
  labels(row1:row2) = j;
end

labels = labels';