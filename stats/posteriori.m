function pp = posteriori(logpri, loglik)
%POSTERIORI Computes the posterior from both log-prior and log-likelihood
%
%   pp = posteriori(logpri, loglik)
%       computes the posterior probabilities of classes given the 
%       logarithm of prior probabilities and the likelihood
%
%       If there are multiple likelihood terms then they can be summed
%       together to form the overall likelihood as
%
%       loglik = loglik1 + loglik2 + ...
%   
%       Suppose there are m classes and n samples, then 
%       - loglik can be either empty (likelihood is missing), 
%         or a matrix of size m x n, where loglik(k, i) gives
%         log p(sample_i | i from model_k).
%       - logpri can be either empty (prior is missing),
%         or a m x 1 column vector (each class has a prior),
%         or a m x n matrix (the prior is sample-dependent).
%
%       In the output, pp will be a matrix of size m x n, where
%       pp(k, i) is p(i from model_k | sample_i).
%
%   Remarks
%       - This function first shifts the value of the log-terms in order
%         to reduce the risk of underflow or overflow.
%

%   History
%   -------
%       - Created by Dahua Lin, on Mar 22, 2009
%       - Modified by Dahua Lin, on Apr 8, 2010
%           - add the support of both sample-dependent and
%             sample-independent prior.
%

%% parse and verify input arguments

assert(isempty(logpri) || isfloat(logpri) && ndims(logpri) == 2, ...
    'posteriori:invalidarg', 'logpri should be a 2D numeric matrix.');

assert(isempty(loglik) || isfloat(loglik) && ndims(loglik) == 2, ...
    'posteriori:invalidarg', 'loglik should be a 2D numeric matrix.');

%% compute

if isempty(logpri)
    L = loglik;
elseif isempty(loglik)
    L = logpri;
else
    if size(logpri, 2) == size(loglik, 2)
        L = loglik + logpri;
    else
        L = bsxfun(@plus, loglik, logpri);
    end
end

L = bsxfun(@minus, L, max(L, [], 1));
pp = exp(L);
pp = bsxfun(@times, pp, 1 ./ sum(pp, 1));

        
