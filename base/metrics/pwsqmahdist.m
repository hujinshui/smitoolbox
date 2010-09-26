function D = pwsqmahdist(X1, X2, A)
% Computes pairwise squared Mahalanobis distances
%
%   D = pwsqmahdist(X1, X2, A);
%       computes pairwise squared Mahalanobis distances between the 
%       columns in X1 and X2.
%
%       The squared Mahalanobis distances are defined as follows
%
%           d(x, y) = (x - y)' * A * (x - y);
%
%       In a d-dimensional space, both X1 and X2 should be respectively 
%       d x n1 and d x n2 matrices, and A should be a d x d symmetric 
%       positive definite matrix.
%
%       In the output, dists will be a n1 x n2 matrix, with 
%           
%           dists(i, j) = dist(X1(:,i), X2(:,j));
%
%       X2 can be input as an empty matrix, which indicates X2 = X1.
%

% Created by Dahua Lin, on Aug 2, 2010
%

%% verify input

if ~(ndims(X1) == 2 && isreal(X1))
    error('pwsqmahdist:invalidarg', 'X1 should be a real matrix.');
end

if ~isempty(X2) && ~(ndims(X1) == 2 && isreal(X2))
    error('pwsqmahdist:invalidarg', 'X2 should be either empty or a real matrix.');
end

d = size(X1, 1);
if ~(ndims(A) == 2 && size(A, 1) == d && size(A, 2) == d)
    error('pwsqmahdist:invalidarg', 'The size of A should be d x d.');
end


%% main

if isempty(X2)
    AX2 = A * X1;
    s1 = sum(X1 .* AX2, 1);
    s2 = s1;
else
    AX2 = A * X2;
    s1 = sum(X1 .* (A * X1), 1);
    s2 = sum(X2 .* AX2, 1);
end
    
D = (-2) * (X1' * AX2);

D = bsxfun(@plus, D, s1');
D = bsxfun(@plus, D, s2);

D(D < 0) = 0;


