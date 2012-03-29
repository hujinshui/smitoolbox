function W = ppca_W(M)
% Get the factor matrix W of a PPCA model
%
%   W = ppca_W(M);
%       computes the factor matrix of the input PPCA model.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% main

W = bsxfun(@times, M.B, M.s);

