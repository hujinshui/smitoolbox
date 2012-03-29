function params = fmm_est_params(pri, gm, X, w, Z)
% Estimate mixture model given sample assignment
%
%   params = fmm_est_params(pri, gm, X, w, Z);
%
%       Estimates component model parameters given sample assignment.
%
%       Input arguments:
%       - gm:       the generative model of observations
%       - pri:      the prior model of parameters
%       - X:        the observations
%       - w:        the weights of observations
%       - Z:        the sample assignment, which can be in either of 
%                   the following forms:
%                   - soft assignment matrix of size K x n
%                   - {K, label_vector}.
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

if isnumeric(Z)
    if isempty(w)
        W = Z.';
    else
        W = bsxfun(@times, Z.', w);
    end 
else
    K = Z{1};
    L = Z{2};    
    W = l2mat(K, L(:), 1, 'sparse');
end

if isempty(pri)
    params = gm.mle(X, W);
else
    cap = gm.capture(X, W);
    params = pri.mapest(cap);
end
