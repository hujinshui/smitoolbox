function V = pwnrmdot(X1, X2)
% Compute pairwise normalized dot products
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
%       - Modified by Dahua Lin, on Aug 2, 2010
%           - simplify error handling
%

%% parse and verify input arguments

if nargin < 2
    X2 = [];
end

if ~(ndims(X1) == 2 && isreal(X1))
    error('pwnrmdot:invalidarg', 'X1 should be a real matrix.');
end

if ~isempty(X2) && ~(ndims(X2) == 2 && isreal(X2))
    error('pwnrmdot:invalidarg', 'X2 should be either empty or a real matrix.');
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


