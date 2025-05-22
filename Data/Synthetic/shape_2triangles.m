function[data,labels] = shape_2triangles(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
rotation = DATAopts.angles;
tau = DATAopts.sigma; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed);

% Re-position each sphere
for j = 1:length(n)

  % Generate noiseless, standardized data
  %thisData = -1+2*rand(n(j),2);
  thisData(:,1) = rand(n(j),1); thisData(:,2) = -1+2*rand(n(j),1);
  thisData = thisData(abs(thisData(:,2)) <= thisData(:,1),:);
  nn(j) = size(thisData,1); 
  thisData(:,d+2:D) = zeros(nn(j),D-d-1); 

  % Generate rotation
  c = cos(rotation(j));
  s = sin(rotation(j));
  G = eye(D);
  G(1:(d+1),1:(d+1)) = [ c  0  s;
                         0  1  0;
                        -s  0  c];

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

  clear thisData
end

labels = labels';