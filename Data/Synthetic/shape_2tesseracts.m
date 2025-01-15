function[data,labels] = shape_2tesseracts(DATAopts)


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
for j = 1:length(n)

  % Generate rotation
  c = cos(rotation(j));
  s = sin(rotation(j));
  G = eye(D);
  G(1:(d+1),1:(d+1)) = [ c 0 0 0 s; 
                         0 1 0 0 0;
                         0 0 1 0 0;
                         0 0 0 1 0;
                        -s 0 0 0 c];
 

  % Generate noiseless, standardized data
  thisData = zeros([n(j) D]);
  thisData(:,1:d) = -0.5 + rand([n(j) d]);

  % Pollute standardized data with noise
  switch noise_type
    case 'normal'
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