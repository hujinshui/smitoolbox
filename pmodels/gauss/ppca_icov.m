function J = ppca_icov(M)
% Get the inverse covariance matrix of a PPCA model
%
%   J = ppca_icov(M);
%       computes the inverse covariance matrix of the input PPCA model.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~is_ppca(M)
    error('ppca_cov:invalidarg', 'M should be a PPCA struct.');
end

%% main

v = 1 / (M.se^2);
a = 1 ./ (1 ./ ((M.s.^2) .* v) + 1);
V = bsxfun(@times, M.B, sqrt(a));
J = v * (eye(M.d) - V * V');

