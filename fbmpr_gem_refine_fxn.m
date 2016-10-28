% FBMPR_GEM_REFINE_FXN   Function for executing the real-valued Fast Bayesian 
% Matching Pursuit algorithm using a repeated greedy search (RGS), accompanied 
% by an iterative Generalized EM (GEM) re-estimation of the signal model
% hyperparameters for the binary (Q = 2) signal model
%
% This function allows one to implement the algorithm described in "Fast 
% Bayesian Matching Pursuit: Model Uncertainty and Parameter Estimation 
% for Sparse Linear Models," P. Schniter, L. C. Potter, and J. Ziniel.  
% FBMPR_GEM_REFINE_FXN will find a high-probability set of sparse 
% solutions "x" to the underdetermined system of linear equations
%                           y = Ax + w,
% where both x and w are real-valued Gaussian.
%
% The refinement stages will progressively refine the given
% hyperparameters using a GEM-based estimation procedure in an attempt to
% converge to hyperparameter choices that minimize NMSE.  This is useful 
% when the hyperparameters are not known with high accuracy.
%
% SYNTAX: [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, d_max, ...
%           e_tot] = fbmpr_gem_refine_fxn(y, A, p1, sig2w, sig2s, mus, ...
%                                        D, stop, E, r_stop)
%
% INPUTS:
%
% y -       Vector of observations of length M
%
% A -       M-by-N measurement matrix 
%
% p1 -      Initial estimate of the prior probability of active taps 
%           (i.e. Pr{s_n = q} = p1/(Q-1), q > 0)
%
% sig2w -   Initial estimate of noise variance.  The authors find that
%           better performance is obtained when the noise is over-estimated
%
% sig2s -   2-by-1 vector of initially estimated sparse coefficient 
%           variances, ([off; on]).  Note that the variances of all Q-1 
%           active Gaussian mixture densities will be equal
%
% mus -     Q-by-1 vector of initially estimated means for the Q-ary 
%           Gaussian mixture densities, ([off; on_1; on_2; ...; on_Q-1]).  
%           Note that the first element of the vector should be zero, to 
%           model sparsity
%
% D -       Maximum number of repeated greedy searches (default 5)
%
% stop -    Threshold for the number of standard deviations below
%           E[nu(s|y)] at which FBMP may terminate upon finding a metric
%           that exceeds it, (default 0)
%
% E -       Maximum number of refinement stages (default 10)
%
% r_stop -  Threshold for early termination of parameter refinement, a
%           number between 0 and 1.  Once the percentage change in
%           parameter estimates for all parameters is below r_stop,
%           FBMP_GEM_REFINE_FXN will terminate, (default 0.02)
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
% d_tot -       Average number of required RGS iterations to find a mixture
%               vector that exceeded E[nu(s,y)] - stop*stddev(nu(s,y)).  If
%               the returned value of d_tot is equal to D, then it is
%               likely that no mixture vectors ever exceeded the threshold
%
% d_max -       RGS iteration at which the maximum metric was found in 
%               the last iteration of the FBMP algorithm, (1 <= d_max <= D)
%
% e_tot -       Number of refinement iterations required, (1 <= e_tot <= E)

%
% Coded by: Justin Ziniel, The Ohio State Univ.
% FBMP algorithm by: Philip Schniter, The Ohio State Univ.
% E-mail: ziniel.1@osu.edu, schniter@ece.osu.edu
% Last change: June 30, 2009
% Change summary: Made minor changes to update this function to FBMPv1.2 
%               function naming conventions.
% FBMP version 1.3 
% Copyright (c) Philip Schniter, Lee C. Potter, and Justin Ziniel, 2008
%

function [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, d_max, e_tot] = ...
                fbmpr_gem_refine_fxn(y, A, p1, sig2w, sig2s, mus, D, stop, E, r_stop)
                            
%% Declarations and error checking

if nargin < 10
    r_stop = 0.02;
    if nargin < 9
        E = 10;
        if nargin < 8
            stop = 0;
            if nargin < 7
                D = 5;
                if nargin < 6
                    error('FMBPR_GEM_REFINE_FXN requires at least 6 arguments')
                end
            end
        end
    end
end

Q = length(mus);
if Q ~= 2
    error(['A GEM-based re-estimation of hyperparameters is only' ...
            ' possible for a binary (Q = 2) signal model'])
end

% Initial parameter estimates
p1_est = p1;
sig2w_est = sig2w;
sig2s_est = sig2s;
mus_est = mus;


%% Begin refinement procedure


for e = 1:E
    
    % Recover signal based on current parameter estimates
    [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot(e), d_max, ...
          hyper_upd] = fbmpr_gem_fxn(y, A, p1_est, sig2w_est, sig2s_est, ...
                                    mus_est, D, stop);
    
    % Transfer estimates from new to old
    p1_old = p1_est;
    sig2w_old = sig2w_est;
    sig2s_old = sig2s_est;
    mus_old = mus_est;
    
    % Refine parameters using newest recovery results
    p1_est = hyper_upd.p1;
    sig2w_est = max([real(hyper_upd.sig2w), 1e-6]);
    sig2s_est = [0; max([real(hyper_upd.sig2s), 1e-2])];
    mus_est = [0; hyper_upd.mu];
    
    % Check for early termination
    percent_change = [abs(p1_est - p1_old)/p1_old; ...
                        abs(sig2w_est - sig2w_old)/sig2w_old; ...
                        abs(sig2s_est(2) - sig2s_old(2))/sig2s_old(2); ...
                        abs(mus_est(2) - mus_old(2))/abs(mus_old(2))];
    if max(percent_change) < r_stop
        break
    end
end

d_tot(d_tot == inf) = D;    % Replace flag by maximum number of iterations
d_tot = mean(d_tot);        % Return the average number of RGS iterations req'd
e_tot = e;                  % Return the total number of refinement stages
