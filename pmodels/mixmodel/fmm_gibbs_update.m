function [S, L] = fmm_gibbs_update(pri, gm, X, w, c0, S, L)
%FMM_GIBBS_UPDATE Gibbs Sampling Update for Finite Mixture Model
%
%   [S, L] = FMM_GIBBS_UPDATE(pri, gm, X, w, S, L);
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

% E-step (re-sampling labels)

Q = ddposterior(S.Pi, L, 'LL');
z = ddsample(Q, 1);
S.z = z;

% M-step (re-estimate component parameters and Pi)

% re-sample parameters

K = S.K;
grps = intgroup(K, z);

params = cell(1, K);
for k = 1 : K
    cap = gm.capture(X, w, grps{k});
    params{k} = pri.pos_sample(cap, 1);
end
params = gm.combine_params(params{:});
S.params = params;

L = gm.loglik(params, X);

% re-sample Pi 

if isempty(w)
    H = intcount(K, z).';
else
    H = aggreg(w, K, z.', 'sum');
end

S.Pi = dird_sample(H + (c0 + 1), 1);



