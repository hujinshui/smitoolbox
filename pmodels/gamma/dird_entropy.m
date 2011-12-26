function v = dird_entropy(alpha, d)
%DIRD_ENTROPY Evaluates the entropy of Dirichlet distribution
%
%   v = dird_entropy(alpha, d);
%
%       Evaluates the entropy of Dirichlet distribution(s).
%       
%       Here, alpha should be a matrix of size d x n, and in the output,
%       the size of v is 1 x n. Particularly, v(i) corresponds to 
%       the parameter given by alpha(:,i).
%
%   v = dird_entropy(alpha, d);
%
%       Evaluates the logarithm of (symmetric) multivariate beta function.
%
%       When alpha is a row vector of size 1 x n, this is equivalent to
%       mvbetaln(repmat(alpha, d, 1)).
%
%       When alpha is a matrix of size d x n, this is equivalent to 
%       mvbetaln(alpha, d).
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% main

if nargin < 2
    d = size(alpha, 1);
end

logB = mvbetaln(alpha, d);

a = alpha - 1;
if size(alpha, 1) == 1
    v = logB + (a * d) .* (psi(alpha * d) - psi(alpha));
else
    salpha = sum(alpha, 1);
    v = logB + sum(a) .* psi(salpha) - sum(a .* psi(alpha));
end

