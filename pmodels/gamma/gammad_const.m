function c = gammad_const(alpha, beta)
% Compute the constant term for gamma distribution logpdf
%
%   c = gammad_const(alpha, beta);
%
%       Evaluates the following constant term for the evaluation of
%       gamma distribution logpdf:
%
%           c = - (alpha log(beta) + gammaln(alpha))
%
%       The sizes of alpha and beta should be compatible in the bsxfun
%       sense.
%
%       If either size(alpha, 1) > 1 or size(beta, 1) > 1, then
%       it treats the distribution as a multi-dimensional distribution.
%
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% main

if isscalar(alpha) || isscalar(beta)    
    c = - (alpha * log(beta) + gammaln(alpha));
else
    c = - bsxfun(@plus, bsxfun(@times, alpha, log(beta)), gammaln(alpha));
end

if size(c, 1) > 1
    c = sum(c, 1);
end

    
    