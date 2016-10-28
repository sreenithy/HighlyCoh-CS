
function [bhat, Sbar, rt] = SWAP(y,X,k,Alg,b,verb)
% Sparse Regression with Correlated Variables
%
% function [bhat, Sbar, rt] = SWAP(y,X,k,Alg,b,verb)
%
% This function solves the following problem:
%   min || y - X * b ||_2^2 s.t. ||b||_0 = k
%
% Inputs:
% y = X * b + noise;
% y <- observations
% X <- measurement matrix
% k <- desired sparsity level in the output (only needed when Alg is
% function)
% Alg <- either a function call or a support set
% b <- the true b vector if known
% verb <- 1 if verbose

% Outputs:
% bhat <- Estimate of b
% Sbar <- Support of bhat
%
% return options
% rt.L = indn;
%
% Example:
%
% p = 1000, n = 200, k = 20;
% X = randn(n,p);
% b = zeros(p,1); b(1:k) = 1;
% y = X * b + randn(n,1);
%
% % Random initialization
% sini = randperm(p); sini = sini(1:k)
% [bhat, Shat] = SWAP(y,X,[],sini,[],0)
%
% % Initialize using OMP
% sini = OMPSupp(y,X,k)
% [bhat, Shat] = SWAP(y,X,[],sini,[],0)
%
% % Initialize using a function call
% % the sparse regression should have the form Alg(y,X,k)
% [bhat, Shat] = SWAP(y,X,k,@OMPSupp,[],0)
%
% author: Divyanshu Vats, dvats@rice.edu
%
% References:
%
% D. Vats and R. G. Baraniuk, "When in Doubt, SWAP: High-Dimensional Sparse
% Recovery from Correlated Measurements," NIPS 2013
%
% D. Vats and R. G. Baraniuk, ""Swapping Variables for High-Dimensional 
% Sparse Regression with Correlated Measurements," Preprint, 2014, arXiv:1312.1706

p = size(X,2);
indn = 1;
ee = .000000000001;

if nargin == 4 || isempty(b)
    b = 0;
end

if ~exist('verb')
    verb = 1;
end

SigmaHat = @(ss,ss1) (X(:,ss)' * X(:,ss1));
omega = X' * y;
LeastSquares = @(S) ( (SigmaHat(S,S) + ee) \ omega(S) );

% Alg is either a support or an algorithm
if strcmp(class(Alg),'function_handle') == 1
    Sbar = Alg(y,X,k);
else
    Sbar = Alg;
end

VV = 1:p;
strue = find(b~=0);

while 1
    %bhat = zeros(size(X,2),1);
    %bhat(Sbar) = LeastSquares(X,Sbar);
    
    SbarC = setdiff(VV,Sbar);
    
    % swap an element from Sbar with an element from SbarC
    Delta = zeros(length(SbarC),length(Sbar));
    %tic;
    for i = 1:length(Sbar)
        Ai = Sbar; Ai(i) = [];
        PiAi = X(:,Ai) * ((SigmaHat(Ai,Ai)) \ X(:,Ai)');
        temp = X(:,Sbar(i)) - PiAi * X(:,Sbar(i));
        ytemp = y' * temp;
        yPi = ytemp * ytemp' / (temp' * temp);
        
        temp = X(:,SbarC) - PiAi * X(:,SbarC);
        ddt = y' * temp;
        ddtt = sum(ddt' .* ddt',2);
        temp_new = sum(temp' .* temp',2);
        
        ddtt = ddtt ./temp_new;
        
        Delta(:,i) = -ddtt + yPi;
    end
    
    [iihat,ihat] = ind2sub(size(Delta),argmin(Delta(:)));
    
    if verb == 1
        display(['NumNon is: ' num2str(length(intersect(strue,Sbar))) ' Delta is: ' num2str(Delta(iihat,ihat)) ]);
    
    end
    
    if Delta(iihat,ihat) > 0
        break;
    end
    
    if length(intersect(strue,Sbar(ihat))) == 1
        display('removed true');
    end
    
    Sbar(ihat) = SbarC(iihat);
    
    indn = indn + 1;
    if indn == p
        break;
    end
    
end

bhat = zeros(p,1);
bhat(Sbar) = LeastSquares(Sbar);
rt.L = indn;

end

function ss = argmin(x)

[t1,t2] = min(x);
ss = t2;

end
