function X = ppca_ft(M, Z)
% Performs forward transform w.r.t. PPCA model
%
%   X = ppca_ft(M, Z);
%       Transforms the latent vectors in Z to the observed space, with 
%       respect to the PPCA model M, as
%
%       x <- W * z.
%
%       Here, Z should be a q x n matrix, and X is a d x n matrix.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~is_ppca(M)
    error('ppca_ft:invalidarg', 'M should be a PPCA struct.');
end

%% main

X = M.B * bsxfun(@times, M.s.', Z);

mu = M.mu;
if ~isequal(mu, 0)
    X = bsxfun(@plus, X, mu);
end

