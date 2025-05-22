function[data,labels] = shape_coneplane(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
tau = DATAopts.sigma; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed); 

labels = [repelem(1,n(1)) repelem(2,n(2))];

% cone
n1 = n(1); X1 = zeros(n(1), D); 
teta=(2*pi)*rand(n1,1);     % Using  spherical coordinates
r=2*rand(n1,1); X1(:,1)=r.*sin(teta); X1(:,2)=r.*cos(teta); X1(:,3)=r;

%plane
n2 = n(2); X2 = zeros(n2, D); 
X2(:,1)=-2+4.*rand(n2,1); X2(:,2)=-2+4.*rand(n2,1); X2(:,3)=1+.2*X2(:,1);

data = [X1; X2]; 

% pollute with noise
% Pollute standardized data with noise
switch noise_type
    case 'gaussian'
      data = data + normrnd(0,tau/sqrt(D-d),size(data));    %noise_stdev/sqrt(D-d)*randn(size(thisData));
    case 'uniform'
      data = data + sqrt(3)*tau/sqrt(D-d)*(2*rand(size(data))-1);
    otherwise
      error("Unknown noise type")
end

end