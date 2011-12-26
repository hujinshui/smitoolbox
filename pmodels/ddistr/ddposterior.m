function Q = ddposterior(pri, lik, op)
%DDPOSTERIOR Compute posterior discrete distribution
%
%   Q = DDPOSTERIOR(pri, lik);
%   Q = DDPOSTERIOR(pri, loglik, 'LL');
%
%       Compute posterior discrete distributions, given prior and
%       likelihood (or log-likelihood).
%
%       Suppose there are K classes and n (observed) samples.
%
%       Input arguments:
%       - pri:      the prior distribution, which can be either of:
%                   - []:   indicates uniform prior
%                   - K x 1 vector: a prior probability vector
%                   - K x n matrix: sample-dependent probability matrix
%
%       - lik:      the likelihood of the samples with respect to 
%                   all K classes, in form of a K x n matrix.
%
%       - loglik:   the log-likelihood of the samples with respect to
%                   all K classes, in form of a K x n matrix.
%
%       Output arguments:
%       - Q:        a K x n matrix, with Q(k, i) being the posterior
%                   probability of the i-th sample for the k-th class.
%

% Created by Dahua Lin, Dec 25, 2011
%

%% verify input arguments

if ~isempty(pri)
    if ~(isfloat(pri) && isreal(pri) && ndims(pri) == 2)
        error('ddposterior:invalidarg', ...
            'pri should be a real vector of matrix.');
    end
end

if ~(isfloat(lik) && isreal(lik) && ndims(lik) == 2)
    error('ddposterior:invalidarg', ...
        'lik or loglik should be a real matrix.');
end


if nargin < 3
    ll = 0;
else
    if ~strcmpi(op, 'll')
        error('ddposterior:invalidarg', ...
            'The third argument can only be ''LL'' if given.');
    end
    ll = 1;
end

[K, n] = size(lik);
if ~isempty(pri)
    if ~(size(pri, 1) == K && (size(pri, 2) == 1 || size(pri, 2) == n))
        error('ddposterior:invalidarg', ...
            'The sizes of pri and lik are inconsistent.');
    end
end

%% main

if K == 1
    Q = ones(1, n);
    return;
end

if ll
    if isempty(pri)
        E = lik;
    else
        if size(pri, 2) == n
            E = log(pri) + lik;
        else
            E = bsxfun(@plus, log(pri), lik);
        end
    end
    Q = nrmexp(E, 1);            
    
else
    if isempty(pri)
        P = lik;
    else
        if size(pri, 2) == n
            P = pri .* lik;
        else
            P = bsxfun(@times, pri, lik);
        end
    end
    Q = bsxfun(@times, P, 1 ./ sum(P, 1));
end

