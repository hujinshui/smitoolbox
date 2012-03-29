function v = gammad_entropy(alpha, beta)
% Compute the entropy for gamma distribution
%
%   v = gammad_entropy(alpha, beta);
%
%       Evaluates the entropy of gamma distribution.
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

u = alpha + gammaln(alpha) + (1 - alpha) .* psi(alpha);

if isscalar(u) || isscalar(beta)
    v = u + log(beta);
else
    v = bsxfun(@plus, u, log(beta));
end

if size(v, 1) > 1
    v = sum(v, 1);
end

