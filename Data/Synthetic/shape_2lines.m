function[data,labels] = shape_2lines(DATAopts)


n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
rotation = DATAopts.angles;
tau = DATAopts.noise_level; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed);

data = zeros([sum(n) D]);
labels = zeros([sum(n) 1]);

% For each line
for j = 1:2

  % Generate rotation
  c = cos(rotation(j));
  s = sin(rotation(j));
  G = eye(D);
  G(1:2,1:2) = [ c  s; 
                -s  c];
 

  % Generate noiseless, standardized data
  thisData = zeros([n(j) D]);
  thisData(:,1) = -0.5 + rand([n(j) 1]);

  % Pollute standardized data with noise
  switch noise_type
    case 'gaussian'
      thisData = thisData + normrnd(0,tau/sqrt(D-d),size(thisData));
    case 'uniform'
      thisData = thisData + sqrt(3)*tau/sqrt(D-d)*(2*rand(size(thisData))-1);
    otherwise
      error("Unknown noise type")
  end

  % Rotate polluted data
  thisData = thisData * G;

  % Place data in right locations in output.
  row1 = 0 + sum(n(1:(j-1))) + 1;
  row2 = sum(n(1:j));
  data(row1:row2,:) = thisData;
  labels(row1:row2) = j;

end