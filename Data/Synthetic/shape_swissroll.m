function[data,labels] = shape_swissroll(DATAopts)

n = DATAopts.number; 
D = DATAopts.ambdim; 
d = DATAopts.intrdim; 
tau = DATAopts.noise_level; 
rng_seed = DATAopts.rngSeed; 
noise_type = DATAopts.noise_type; 

rng(rng_seed); 

% The roll 
lb=pi; ub = 4*pi; ub2 = 40; n1 = n(1); %74.7094;
counter = 0; Z = zeros(n1,1);
while counter<n1
    possible_point = [ub*rand ub*rand];
    if possible_point(2) < possible_point(1) && possible_point(1) > lb
        Z(counter+1)=possible_point(1);
        counter = counter+1;
    end
end

X1 = zeros(n1,3);
for i=1:n1
   X1(i,:) = [3*Z(i)*cos(Z(i))/ub2 rand 3*Z(i)*sin(Z(i))/ub2]; 
end
X1(:,d+2:D) = zeros(n1,D-d-1); 

% the plane
n2=n(2);X2(:,1)=-1+2*rand(n2,1); X2(:,2)=rand(n2,1);X2(:,d+1:D) = zeros(n2,D-d); 

data = [X1; X2]; %data(:,d+2:D) = zeros(sum(n),D-d-1);
% Pollute standardized data with noise
switch noise_type
    case 'gaussian'
      data = data + normrnd(0,tau/sqrt(D-d),size(data));    %noise_stdev/sqrt(D-d)*randn(size(thisData));
    case 'uniform'
      data = data + sqrt(3)*tau/sqrt(D-d)*(2*rand(size(data))-1);
    otherwise
      error("Unknown noise type")
end

labels = [repelem(1,n(1)) repelem(2, n(2))]';
end