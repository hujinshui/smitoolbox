function [S, L] = fmm_em_update(pri, gm, X, w, c0, S, L)
%FMM_EM_UPDATE Expectation-Maximization Update for Finite Mixture Model
%
%   [S, L] = FMM_EM_UPDATE(pri, gm, X, w, S, L);
%
%       Updates finite mixture model solution as an E-M iteration.
%
%       Input arguments:
%       - pri:      The prior object
%       - gm:       The generative model object
%       - X:        The samples
%       - w:        The sample weights
%       - c0:       The prior count of components
%
%       The arguments to be updated:
%       - S:        The finite-mixture model solution.
%       - L:        The log-likelihood matrix.
%
%   Remarks
%   -------
%       - This function is to help the implementation of various mixture
%         model estimation, and is not supposed to be directly called by
%         end users.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% main

% E-step (re-estimate marginal probabilities)

Q = ddposterior(S.Pi, L, 'LL');
S.Q = Q;

% M-step (re-estimate component parameters and Pi)

% estimate component parameters

params = fmm_est_params(pri, gm, X, w, Q);

S.params = params;
L = gm.loglik(params, X);

% estimate Pi

S.Pi = ddestimate(Q, w, c0);



