function X = ppca_sample(M, n)
% Draws n samples from a PPCA model
%
%   X = ppca_sample(M, n);
%       draws n samples from a PPCA model M.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~is_ppca(M)
    error('ppca_sample:invalidarg', 'M should be a ppca struct.');
end

if nargin < 2
    n = 1;
else
    if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 1)
        error('ppca_sample:invalidarg', 'n should be a positive integer.');
    end
end

%% main

Z = randn(M.q, n);
X = ppca_ft(M, Z) + randn(M.d, n) * M.se;

