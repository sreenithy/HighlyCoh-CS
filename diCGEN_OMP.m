clear all; 
close all; 
display('===========================================')

noise    = 5;   % noise level ||z||_2/||\sum_{j=1}^s x_j*e^{-2*\pi*i*k*supp(j)}||_2 = noise/100
range    = 10;  %dyanmic range = xmax/xmin
fc   = 25;
RL   = 1/(2*fc);   %Rayleigh length
sepn1= 2;
sepn2= 3;
sepn = sepn1;
sep1 = sepn1*RL;
sep2 = sepn2*RL;
supp = sep1;
while supp(end) < 1-2*sep2
    move = sep1+rand*(sep2-sep1);
    supp = [supp ; supp(end)+move];
end
sparsity = length(supp);

k         = -fc:1:fc;
n         = 2*fc+1;
amp  = rand(sparsity,1);
amp  = (amp-min(amp))/(max(amp)-min(amp));      %[0,1]
amp  = 1 + amp*(range-1);
sig  = (binornd(1,1/2,sparsity,1)-0.5)*2;
x    = amp.*sig;%.*exp(2*pi*1i*rand(sparsity,1));
xabs = abs(x);
F = exp(-1i*2*pi*k'*supp'); % Fourier matrix
y = F*x;
sigma    = noise/100*norm(y)/sqrt(2*n);
z        = normrnd(0,sigma,n,1) + 1i*normrnd(0,sigma,n,1);
y        = y + z;
if noise > 0                       % ||z|| <= ep
    ep = sqrt(2*n)*sigma;
else
    ep = 10^(-12);
end
fprintf('noise to signal ratio = %6.2fp \n',100*norm(z)/norm(y))
N = 2*fc+1;            % the number of rows of sensing matrix A
F = 10;                % refinement factor / super-resolution factor
M = N * F;             % the number of columns of sensing matrix A
%% construct the sensing matrix A
A = exp(-1i*2*pi*(-fc:fc)'*(0:M-1)/M);  
size(A)
[N,M] = size(A);

cohMatrix = zeros(M,M);
for i = 1 : M
    for j = i+1 : M
        cohMatrix(i,j) = A(:,i)'*A(:,j)/(norm(A(:,i)) * norm(A(:,j)));
    end
end

minCoh = min(min(cohMatrix));
maxCoh = max(max(cohMatrix))
cohmatrix = cohMatrix + cohMatrix' + eye(M);
