function dists = eucdist(X1, X2, var3, var4)
%Computes the Euclidean distances between corresponding vectors.
%
%   dists = eucdist(X1, X2);
%       computes the Euclidean distances between corresponding vectors in
%       X1 and X2.
%
%       The Euclidean distance, also known as L2 distance, is defined as
%
%           d(x, y) = sqrt( sum_i ( x(i) - y(i) )^2 )
%
%       X1 and X2 in the input should be of the same size. Suppose they
%       are both d x n matrices, then dists will be a 1 x n vector, with
%
%           dists(i) = d(X1(:, i), X2(:, i));
%
%   dists = eucdist(X1, X2, 'square');
%       computes the square of Euclidean distances.
%
%       The computation is basically the same, except that the sum of
%       difference squares rather than its square root is returned.
%
%   dists = eucdist(X1, X2, w);
%       computes the weighted Euclidean distances, which is defined as
%
%           d(x, y) = sqrt( sum_i w(i) * ( x(i) - y(i) )^2 )
%
%       w should be an d x 1 column vector specifying the weights of
%       different components, all weights should be non-negative.
%
%   dists = eucdist(X1, X2, w, 'square');
%       computes the weighted squared Euclidean distances.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'eucdist:invalidarg', ...
    'X1 and X2 should be both numeric matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'eucdist:invalidsize', ...
    'X1 and X2 should be of the same size.');

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
    assert(isnumeric(w) && ndims(w) == 2 && size(w,1) == d && size(w,2) == 1, ...
        'eucdist:invalidarg', ...
        'w should be a d x 1 numeric vector.');
    
    assert(all(w >= 0), ...
        'eucdist:negativeweight', ...
        'all weights should be non-negative.');
end


%% main

D = X1 - X2;

if ~weighted
    D = sum(D .* D, 1);
else
    D = sum(bsxfun(@times, D .* D, w), 1);
end

if sqr
    dists = D;
else
    dists = sqrt(D);
end


