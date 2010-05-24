function V = pwnrmdot(X1, X2)
%PWNRMDOT Computes pairwise normalized dot products
%
%   V = pwnrmdot(X1, X2);
%       computes pairwise dot products between column vectors in X1 and X2.
%
%       Normalized dot-product is defined as the dot-product between
%       vectors after L2-normalization, as
%       
%           normalized dot-product(x, y) = x' * y / (||x|| * ||y||)
%
%       Suppose X1 and X2 are respectively d x n1 and d x n2 matrices, then
%       V will be a n1 x n2 matrix, with V(i, j) being the normalized dot
%       product between X1(:, i) and X2(:, j).
%
%   V = pwnrmdot(X);
%       computes pairwise dot products between columns in X, which is
%       functionally equivalent to pwnrmdot(X, X).
%
%       However, it takes advantages of the fact that X1 = X2 = X to make
%       slightly faster implementation.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

if nargin < 2
    assert(ndims(X1) == 2, 'pwnrmdot:invalidarg', ...
        'X1 should be a matrix.');
    
    X2 = [];
    
elseif ~isempty(X2)    
    assert(ndims(X1) == 2 && ndims(X2) == 2 && size(X1,1) == size(X2,1), ...
        'pwnrmdot:invalidarg', ...
        'X1 and X2 should be both matrices with the same number of rows.');
    
end

%% main

if isempty(X2)
   
    n = size(X1, 2);
    
    V = X1' * X1;
    s = V(1 + (n+1) *(0:n-1));  % take the diagonal elements
    
    q = 1 ./ sqrt(s);
    V = bsxfun(@times, V, q);
    V = bsxfun(@times, V, q');        
    
else
    
    s1 = sum(X1 .* X1, 1);
    s2 = sum(X2 .* X2, 1);
    
    V = X1' * X2;
    V = bsxfun(@times, V, 1 ./ sqrt(s1)');
    V = bsxfun(@times, V, 1 ./ sqrt(s2));
        
end

V(V > 1) = 1;
V(V < -1) = -1;


