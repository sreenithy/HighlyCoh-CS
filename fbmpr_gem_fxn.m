% FBMPR_GEM_FXN   Function for executing the real-valued Fast Bayesian Matching 
% Pursuit algorithm using a repeated greedy search (RGS) followed by a 
% generalized EM (GEM) re-estimation of hyperparameters for a binary (Q = 2) 
% signal model
%
% This function allows one to implement the algorithm described in "Fast 
% Bayesian Matching Pursuit: Model Uncertainty and Parameter Estimation 
% for Sparse Linear Models," P. Schniter, L. C. Potter, and J. Ziniel.  
% FBMPR_GEM_FXN will find a high-probability set of sparse solutions "x"
% to the underdetermined system of linear equations
%                           y = Ax + w,
% where both x and w are real-valued Gaussian.
%
% SYNTAX: [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, ...
%          d_max, hyper_upd] = fbmpr_gem_fxn(y, A, p1, sig2w, sig2s, ...
%                                           mus, D, stop)
%
% INPUTS:
%
% y -       Vector of observations of length M
%
% A -       M-by-N measurement matrix 
%
% p1 -      Prior probability of active taps (i.e. Pr{s_n = q} = p1/(Q-1),
%           q > 0)
%
% sig2w -   Noise variance
%
% sig2s -   2-by-1 vector of sparse coefficient variances,  ([off; on]).  
%           Note that the variances of all Q-1 active Gaussian mixture
%           densities will be equal
%
% mus -     Q-by-1 vector of means for the Q-ary Gaussian mixture 
%           densities, ([off; on_1; on_2; ...; on_Q-1]).  Note that the 
%           first element of the vector should be zero, to model sparsity.
%
% D -       Maximum number of repeated greedy searches (default 5)
%
% stop -    Threshold for the number of standard deviations below
%           E[nu(s|y)] at which FBMP may terminate upon finding a metric
%           that exceeds it, (default 0)
%
%
% OUTPUTS:
%
% xmmse -       Approximate sparse MMSE estimate of x, given y
%
% xmmse_star -	Significant conditional MMSE estimates of x, given s and y
%
% psy_star -    Corresponding significant conditional mixture vector
%               posterior probabilities, p(s|y)
%
% nu_star -     Corresponding metric values, nu(s,y)
%
% T_star -      Indices of active taps for the conditional MMSE estimates 
%               of x returned in xmmse_star
%
% d_tot -       Number of required RGS iterations to find a mixture vector
%               that exceeded E[nu(s,y)] - stop*stddev(nu(s,y)).  If no 
%               such vector is found, 'inf' is returned
%
% d_max -       RGS iteration at which the maximum metric was found, 
%               (1 <= d_max <= d_tot)
%
% hyper_upd -   The updated values of the hyperparameters for the binary (Q
%               = 2) signal model, returned in a structure with the
%               following elements:
%                   .p1 - Updated estimate of p1, Pr{s_n ~= 0}
%                   .sig2w - Updated estimate of the noise variance
%                   .sig2s - Updated estimate of active coefficient
%                            variance
%                   .mu - Updated estimate of active coefficient mean

%
% Coded by: Philip Schniter, The Ohio State Univ.
% Modified by: Justin Ziniel, The Ohio State Univ.
% E-mail: schniter@ece.osu.edu, ziniel.1@osu.edu
% Last change: June 30, 2009
% Change summary: Corrected beta definition and modified code to allow for
%                 complex-valued A matrices and observations (1/17/09)
%               Allow A with non-unit-norm columns and update to 
%               FBMPv1.2 function naming conventions (6/30/09)
% FBMP version 1.3 
% Copyright (c) Philip Schniter, Lee C. Potter, and Justin Ziniel, 2008
%

function [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, d_max, ...
          hyper_upd] = fbmpr_gem_fxn(y, A, p1, sig2w, sig2s, mus, D, stop)


%% Fill in missing variables and check input for errors

% Set defaults
if (nargin < 8)
    stop = 0;           % Stop at zero std. devs. below the mean of nu(s,y)
    if (nargin < 7)
        D = 5;          % Maximum number of RGS iterations
        if (nargin < 6)
            error('FBMPR_GEM_FXN requires at least 6 arguments')
        end
    end
end

[M, N] = size(A);

if (length(y) ~= M)
    error('y and A have a differing number of rows')
elseif (p1 < 0 || p1 > 1)
    error('Invalid sparseness prior, p1')
elseif (length(sig2s) ~= 2)
    error('Vector of mixture density variances must be of length 2')
elseif (min(sig2s) < 0) || (sig2w < 0)
    error('Variances must be non-negative')
elseif (D < 1)
    error('D < 1')
elseif (length(mus) > 2)
    error(['The generalized EM update of hyperparameters is only '...
          'compatible with a binary (Q = 2) signal model'])
end

Q = length(mus) - 1;        % # of Gaussian mixture densities minus one
if Q ~= 1
    error(['A GEM re-estimation of hyperparameters is only possible'...
            ' for a binary (Q = 2) signal model'])
end
ps = [1 - p1; ones(Q,1)*p1/Q];      % All active Gaussians are equiprobable
sig2s = [sig2s; sig2s(2)*ones(Q-1,1)];  % Equal variance for all active Gaussians

% 1st and 2nd moments of mixture selection metric prior distribution
a2 = trace(A'*A)/N;			% average 2-norm of the columns of A
nu_true_mean = -M/2 - M/2*log(sig2w) - p1*N/2*log(a2*sig2s(2)/sig2w+1) - ...
                    M/2*log(2*pi) + N*log(ps(1)) + p1*N*log(ps(2)/ps(1));
nu_true_stdv = sqrt(M/2 + N*p1*(1-p1)*(log(ps(2)/ps(1)) - ...
                    log(a2*sig2s(2)/sig2w + 1)/2)^2);
nu_stop = nu_true_mean - stop*nu_true_stdv;


% Default algorithmic parameters
psy_thresh = 1e-4;		% significant-posterior threshold 
P = min(M, 1 + ceil(N*p1 + erfcinv(1e-2)*sqrt(2*N*p1*(1 - p1))));  % search length
H = 10;     % Total number of "M" steps to perform during the GEM hyperparameter
            % estimation procedure
gem_keep = 10;   % How many mixture vectors of S_star to use in the GEM
                % update procedure


%% Repeated greedy search

% Allocate variables for storage
T = cell(P,D);          % indices of active taps
sT = cell(P,D);         % active mixing params
nu = -inf*ones(P,D);	% metrics
xmmse = cell(P,D);		% mmse conditioned on active taps
d_tot = inf;            % flag maximum number of RGS iterations at infinity

% Initialize (root node)
nu_root = -norm(y)^2/2/sig2w - M*log(2*pi)/2 - M*log(sig2w)/2 + N*log(ps(1));
Bxt_root = A/sig2w;
betaxt_root = abs(sig2s(2)*(1 + sig2s(2)*sum(conj(A).*Bxt_root)).^(-1));
nuxt_root = zeros(1,Q*N);	% Q = # non-zero means
for q = 1:Q
    nuxt_root([1:N]+(q-1)*N) = nu_root + log(betaxt_root/sig2s(2))/2 ...
            + 0.5*betaxt_root.*abs( y'*Bxt_root + mus(q + 1)/sig2s(2)).^2 ...
            - 0.5*abs(mus(q + 1))^2/sig2s(2) + log(ps(2)/ps(1));
end

% Descend one branch at a time
for d = 1:D,
    nuxt = nuxt_root;
    z = y;
    Bxt = Bxt_root;
    betaxt = betaxt_root;
    for p = 1:P
        [nustar,nqstar] = max(nuxt);                % find best extension
        while sum(abs(nustar-nu(p,1:d-1)) < 1e-8)   % if same as explored node...
            nuxt(nqstar) = -inf;                    % ... mark extension as redundant
            [nustar, nqstar] = max(nuxt);           % ... and find next best extension
        end
        qstar = floor((nqstar - 1)/N) + 1;          % mean index of best extension
        nstar = mod(nqstar - 1, N) + 1;             % coef index of best extension
        nu(p,d) = nustar;                           % replace worst explored node...
        if (p > 1)
            T{p,d} = [T{p-1,d}, nstar];
            sT{p,d} = [sT{p-1,d}, 1 + qstar];
        else
            T{p,d} = nstar;
            sT{p,d} = 1 + qstar;
        end
        z = z - A(:,nstar)*mus(qstar+1);
        Bxt = Bxt - Bxt(:,nstar)*betaxt(nstar)*( Bxt(:,nstar)'*A );
        xmmse{p,d} = zeros(N,1);
        xmmse{p,d}(T{p,d}) = mus(sT{p,d}) + sig2s(2)*Bxt(:,T{p,d})'*z;
        betaxt = abs(sig2s(2)*(1 + sig2s(2)*sum(conj(A).*Bxt)).^(-1));
        for q = 1:Q                                 % Q = # non-zero means
            nuxt([1:N]+(q-1)*N) = nu(p,d) + log(betaxt/sig2s(2))/2 ...
                    + 0.5*betaxt.*abs( z'*Bxt+mus(q+1)/sig2s(2) ).^2 ...
                    - 0.5*abs(mus(q+1))^2/sig2s(2)  + log(ps(2)/ps(1));
            % can't activate an already activated coefficient!
            nuxt(T{p,d}+(q-1)*N) = -inf*ones(size(T{p,d}));
        end
    end

    if (max(nu(:,d)) > nu_stop)     % A mixture vector has exceeded the threshold
        d_tot = d;
        break
    end
end
nu = nu(:,1:d);


%% Calculate stuff for statistically significant hypotheses

[dum,indx] = sort(nu(:), 'descend');
d_max = ceil(indx(1)/P);
nu_max = nu(indx(1));
num = sum( nu(:) > nu_max+log(psy_thresh) );
nu_star = nu(indx(1:num));
psy_star = exp(nu_star - nu_max)/sum(exp(nu_star - nu_max));    % posteriors
T_star = cell(1,num);
xmmse_star = cell(1,num);           % mmse est cond on mixing params
for k = 1:num
    T_star{k} = T{indx(k)};
    xmmse_star{k} = xmmse{indx(k)};
end

xmmse = [xmmse_star{:}]*psy_star;  % approximate mmse estimate


%% GEM-based re-estimate of binary signal model hyperparameters

alpha = sig2w/sig2s(2);     % Initial value of NSR
gem_keep = min([num, gem_keep]);    % Make sure there are at least gem_keep terms
psy_min = psy_star(1:gem_keep);     % Truncated terms to keep
psy_min = (1/sum(psy_min))*psy_min;	% Normalized
for h = 1:1:H       % Number of coordinate ascent steps to perform on sig2s
    % MMSE-based hyperparameter re-estimation
    if h == 1       % Only perform one "M" step
        hyper_upd.p1 = 0;                   % Initialize hyperparameter updates
        mu_numerator = 0;
        mu_denominator = 0;
        mu_multiplier = zeros(M,1);
        hyper_upd.sig2w = 0;
        for k = 1:gem_keep
            s = zeros(N,1);
            s(T_star{k}) = 1;
            A_s = A(:,T_star{k});
            cov_inv = (1/alpha)*(eye(M) - (A_s/(alpha*eye(sum(s)) + ...
                        A_s'*A_s))*A_s');   % MIL inverse
            hyper_upd.p1 = hyper_upd.p1 + psy_min(k)*sum(s);
            mu_multiplier = psy_min(k)*s'*A'*cov_inv;
            mu_numerator = mu_numerator + mu_multiplier*y;
            mu_denominator = mu_denominator + mu_multiplier*A*s;
        end
        hyper_upd.p1 = (1/N)*hyper_upd.p1;              % p1 update
        hyper_upd.mu = mu_numerator/mu_denominator;     % mu update
    end

    hyper_upd.sig2s = 0;                % Initialize hyperparameter update
    for k = 1:gem_keep
        s = zeros(N,1);
        s(T_star{k}) = 1;
        A_s = A(:,T_star{k});
        cov_inv = (1/alpha)*(eye(M) - (A_s/(alpha*eye(sum(s)) + ...
                    A_s'*A_s))*A_s');   % MIL inverse
        hyper_upd.sig2s = hyper_upd.sig2s + ...
                psy_min(k)*((y - hyper_upd.mu*A*s)'*cov_inv)...
                *(y - hyper_upd.mu*A*s);
        if h == 1       % Only perform one "M" step
            cov_inv = (eye(M) - (A_s/(A_s'*A_s))*A_s');
            hyper_upd.sig2w = hyper_upd.sig2w + (1/(M - sum(s)))*psy_min(k)*...
                    ((y - hyper_upd.mu*A*s)'*cov_inv)*...
                    (y - hyper_upd.mu*A*s);
        end
    end
    % Multiply GEM-based re-estimate by needed constants
    hyper_upd.sig2s = (1/M)*hyper_upd.sig2s;
    alpha = hyper_upd.sig2w/hyper_upd.sig2s;     % Updated value of NSR
end
