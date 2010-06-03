function Y = mmvmult(A, X)
%MMVMULT Performs multiple matrix-vector multiplication
%
%   Y = mmvmult(A, X);
%       performs multiple matrix-vector multiplication.
%       Here, A should be an array of size p x q x n, and X of size q x n,
%       then Y is of size p x n, such that
%
%           Y(:,i) = A(:,:,i) * X(:,i).
%

% Created by Dahua Lin, on Apr 4, 2010
%

%% verify input arguments

assert(isfloat(A) && ndims(A) <= 3, 'mmvmult:invalidarg', ...
    'A should be a numeric array with ndims(A) <= 3');

[p, q, n] = size(A);

assert(isfloat(X) && isequal(size(X), [q n]), 'mmvmult:invalidarg', ...
    'X should be a numeric matrix of size q x n.');

%% main

if p == 1
    Y = sum(reshape(A, q, n) .* X, 1);
    
elseif q == 1
    Y = bsxfun(@times, reshape(A, p, n), X);
        
elseif p < n && p < q  % small output dimension
    Y = zeros(p, n, class(A(1) * X(1)));    
    for k = 1 : p
        Y(k, :) = sum(reshape(A(k,:,:), q, n) .* X, 1);
    end
    
elseif q < n && q < p  % small input dimension
    Y = bsxfun(@times, reshape(A(:,1,:), p, n), X(1, :));
    for k = 2 : q
        Y = Y + bsxfun(@times, reshape(A(:,k,:), p, n), X(k, :));
    end    
    
else  % large output dimension 
    Y = zeros(p, n, class(A(1) * X(1))); 
    for i = 1 : n
        Y(:,i) = A(:,:,i) * X(:,i);
    end    
end
