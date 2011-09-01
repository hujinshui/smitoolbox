function v = gamma_entropy(alpha)
% Computes the entropy of gamma distribution(s)
%
%   v = gamma_entropy(alpha);
%
%       Computes the entropy of a Gamma distribution with shape parameter
%       alpha and unit scale.
%
%       The output array has the same size as alpha. 
%

% Created by Dahua Lin, on Sep 1, 2010
%

%% verify input

if ~(isfloat(alpha))
    error('gamma_entropy:invalidarg', 'alpha should be a numeric array.');
end

%% main

v = alpha + gammaln(alpha) + (1 - alpha) .* psi(alpha);
