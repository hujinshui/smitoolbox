function v = nrmdot(X1, X2)
%NRMDOT Computes the normalized dot-product between corresponding vectors
%
%   v = nrmdot(X1, X2);
%       computes the normalized dot-product between corresponding column 
%       vectors in X1 and X2.
%
%       Normalized dot-product is defined as the dot-product between
%       vectors after L2-normalization, as
%       
%           normalized dot-product(x, y) = x' * y / (||x|| * ||y||)
%
%       X1 and X2 should be matrices of the same size. Let their size be
%       d x n, then v will be a 1 x n vector with v(i) being the normalized
%       dot product between X1(:, i), and X2(:, i);
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Jul 22, 2010
%           - simplify error handling
%

%% parse and verify input arguments

if ~(isfloat(X1) && isfloat(X2) && ndims(X1) == 2 && ndims(X2) == 2)
    error('nrmdot:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
end


%% main

s1 = sum(X1 .* X1, 1);
s2 = sum(X2 .* X2, 1);

v = sum(X1 .* X2, 1) ./ sqrt(s1 .* s2);

