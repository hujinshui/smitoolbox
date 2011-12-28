function Z = ppca_bt(M, X)
% Performs backward transform w.r.t. PPCA model
%
%   X = ppca_ft(M, Z);
%       Transforms the observed vectors in X to the latent space, with 
%       respect to the PPCA model M, as
%
%       Here, X should be a d x n matrix, and Z is a q x n matrix.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~is_ppca(M)
    error('ppca_bt:invalidarg', 'M should be a PPCA struct.');
end

%% main

mu = M.mu;
if ~isequal(mu, 0)
    X = bsxfun(@minus, X, mu);
end

Z = bsxfun(@times, M.B' * X, 1 ./ M.s.');

