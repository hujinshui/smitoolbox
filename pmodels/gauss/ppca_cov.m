function C = ppca_cov(M)
% Get the covariance matrix of a PPCA model
%
%   C = ppca_cov(M);
%       computes the equivalent covariance matrix of the input PPCA model.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~is_ppca(M)
    error('ppca_cov:invalidarg', 'M should be a PPCA struct.');
end

%% main

W = ppca_W(M);
C = adddiag(W * W', M.se^2);

