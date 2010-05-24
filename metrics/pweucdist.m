function dists = pweucdist(X1, X2, var3, var4)
%PWEUCDIST Computes pairwise Euclidean distances
%
%   dists = pweucdist(X1, X2);
%       computes the Euclidean distances between the columns in X1 and X2.
%
%       The Euclidean distance, also known as L2 distance, is defined as
%
%           d(x, y) = sqrt( sum_i ( x(i) - y(i) )^2 )
%
%       Suppose X1 and X2 are respectively d x n1 matrix and d x n2 matrix,
%       then dists is a matrix of size n1 x n2, such that dists(i, j) is 
%       the Euclidean distance between X1(:, i) and X2(:, j).
%
%   dists = pweucdist(X1, X2, 'square');
%       computes the square of pairwise Euclidean distances.
%
%       The computation is basically the same, except that the sum of
%       difference squares rather than its square root is returned.
%
%   dists = pweucdist(X1, X2, w);
%       computes the pairwise weighted Euclidean distances, which is defined as
%
%           d(x, y) = sqrt( sum_i w(i) * ( x(i) - y(i) )^2 )
%
%       w should be an d x 1 column vector specifying the weights of
%       different components, all weights should be non-negative.
%
%   dists = pweucdist(X1, X2, w, 'square');
%       computes the square of pairwise weighted Euclidean distances.
%       
%
%   The following simplified syntax are also supported.
%
%   dists = pweucdist(X);
%   dists = pweucdist(X, []);
%   dists = pweucdist(X, [], 'square');
%   dists = pweucdist(X, [], w);
%   dists = pweucdist(X, [], w, 'square');
%       computes the pairwise Euclidean distances (or their squares) for
%       the columns in X.
%
%       By explicitly telling the function that X2 = X1 via an empty matrix
%       input at X2, a faster implementation can be used.
%       

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Jun 3, 2008
%           - implement a faster algorithm for the case that X2 = [].
%


%% parse and verify input arguments

if nargin < 2 
    X2 = [];    
end

if isempty(X2)
    assert(isnumeric(X1) && ndims(X1) == 2, ...
        'pweucdist:invalidarg', ...
        'X1 should be a numeric matrix.');
else
    assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
        'pweucdist:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
    
    assert(size(X1, 1) == size(X2, 1), ...
        'pweucdist:invalidsize', ...
        'Columns in X1 and X2 differ in length.');
end


if nargin < 3
    sqr = false;
    weighted = false;
    
elseif nargin == 3
    if isnumeric(var3)
        sqr = false;
        weighted = true;
        w = var3;
        
    elseif strcmp(var3, 'square')
        sqr = true;
        weighted = false;
        
    else
        error('eucdist:invalidarg', ...
            'The 3rd input argument is invalid.');
        
    end
    
else  % nargin == 4
    if isnumeric(var3) && strcmp(var4, 'square')
        sqr = true;
        weighted = true;
        w = var3;
    else
        error('eucdist:invalidarg', ...
            'The 3rd or 4th input argument is invalid.');
    end
end

if weighted
    d = size(X1, 1);
    
    assert(isnumeric(w) && ndims(w) == 2 && size(w,1) == d && size(w,2) == 1, ...
        'eucdist:invalidarg', ...
        'w should be a d x 1 numeric vector.');
    
    assert(all(w >= 0), ...
        'eucdist:negativeweight', ...
        'all weights should be non-negative.');
end


%% main

if isempty(X2)
   
    if weighted
        X1 = bsxfun(@times, X1, sqrt(w));        
    end
        
    n = size(X1, 2);
    D = X1' * X1;
    
    % take the diagonal elements, which equals sum(X1 .* X1, 1);
    sx = D(1 + (n+1) * (0:n-1));  
    
    D = bsxfun(@plus, (-2) * D, sx);
    D = bsxfun(@plus, D, sx');
    
else
    
    if weighted
        sw = sqrt(w);
        X1 = bsxfun(@times, X1, sw);
        X2 = bsxfun(@times, X2, sw);
    end
    
    D = (-2) * (X1' * X2);            
    D = bsxfun(@plus, D, sum(X1 .* X1, 1)');
    D = bsxfun(@plus, D, sum(X2 .* X2, 1));
        
end


% enforce non-negativeness (rounding error may lead to negative values)
D(D < 0) = 0;

if sqr
    dists = D;
else
    dists = sqrt(D);
end

