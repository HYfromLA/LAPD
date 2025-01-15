function[data,labels] = shape_dollarsign(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
tau = DATAopts.noise_level; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed);

% For each manifold

% The S shape. 
n1 = n(1); 
X1 = randn(n1, d+1); X1 = X1 ./ vecnorm(X1,2,2);
X11 = X1(1:floor(0.5*n1),:); X11(:,2) = X11(:,2)+1; X11(X11(:,1)>0 & X11(:,2)<1,:)=[];
X12 = X1(floor(0.5*n1+1):n1,:); X12(:,2) = X12(:,2)-1; X12(X12(:,1)<0 & X12(:,2)>-1,:)=[];

% The bar. 
X1 = [X11;X12]; X2 =[repelem(0,n(2))' 5*rand(n(2), 1)-2.5]; 

data = [X1; X2]; data(:,d+2:D) = zeros(size(data,1),D-d-1); 
 
% Pollute standardized data with noise
switch noise_type
    case 'gaussian'
      data = data + normrnd(0,tau/sqrt(D-d),size(data));    %noise_stdev/sqrt(D-d)*randn(size(thisData));
    case 'uniform'
      data = data + sqrt(3)*tau/sqrt(D-d)*(2*rand(size(data))-1);
    otherwise
      error("Unknown noise type")
end

labels = [repelem(1, size(X1,1)) repelem(2, n(2))]';

end