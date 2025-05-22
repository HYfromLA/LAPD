function[data,labels] = shape_rosecircles(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
tau = DATAopts.sigma; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed);

% For each manifold

% The rose. 
n1 = n(1);
theta=(2*pi)*randn(n1,1);    
x1=2*cos(theta/.5).*sin(theta);    
y1=2*cos(theta/.5).*cos(theta);   
X1=[x1 y1]; %X1(:,d+2:D) = zeros(n(1),D-d-1); 

% The two circles. 
X2 = randn(n(2), 2); X2 = X2 ./ vecnorm(X2,2,2);
X3 = randn(n(3), 2); X3 = X3 ./ vecnorm(X3,2,2)+1;

data = [X1; X2; X3]; data(:,d+2:D) = zeros(sum(n),D-d-1); 
 
% Pollute standardized data with noise
switch noise_type
    case 'gaussian'
      data = data + normrnd(0,tau/sqrt(D-d),size(data));    %noise_stdev/sqrt(D-d)*randn(size(thisData));
    case 'uniform'
      data = data + sqrt(3)*tau/sqrt(D-d)*(2*rand(size(data))-1);
    otherwise
      error("Unknown noise type")
end

labels = [repelem(1,n(1)) repelem(2, n(2)) repelem(3, n(3))]';

end