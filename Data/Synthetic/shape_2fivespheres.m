function[data,labels] = shape_2fivespheres(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
tau = DATAopts.sigma; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed);

data = zeros([sum(n) D]);
labels = zeros([sum(n) 1]);

% Re-position each sphere
for j = 1:length(n)

  % Generate noiseless, standardized data
  thisData = randn(n(j), d+1); 
  thisData = thisData ./ vecnorm(thisData,2,2);
  if j>1
    thisData(:,2:d) = thisData(:,2:d) + 1.5/sqrt(d-1);
  end
  thisData(:,d+2:D) = zeros(n(j),D-d-1); 

  % Pollute standardized data with noise
  switch noise_type
    case 'gaussian'
      thisData = thisData + normrnd(0,tau/sqrt(D-d),size(thisData));    %noise_stdev/sqrt(D-d)*randn(size(thisData));
    case 'uniform'
      thisData = thisData + sqrt(3)*tau/sqrt(D-d)*(2*rand(size(thisData))-1);
    otherwise
      error("Unknown noise type")
  end

  % Place data in right locations in output.
  row1 = 0 + sum(n(1:(j-1))) + 1;
  row2 = sum(n(1:j));
  data(row1:row2,:) = thisData;
  labels(row1:row2) = j;
end