function dists = cityblkdist(X1, X2, w)
%CITYBLKDIST Computes the city-block distances between corresponding vectors
%
%   dists = cityblkdist(X1, X2);
%       computes the city-block distances between corresponding vectors in
%       X1 and X2.
%
%       The city-block distance, also known as L1 distance, rectilinear 
%       distance, or Manhattan distance, is defined as
%
%           d(x, y) =  sum_i | x(i) - y(i) |
%
%       X1 and X2 in the input should be of the same size. Suppose they
%       are both d x n matrices, then dists will be a 1 x n vector, with
%
%           dists(i) = d(X1(:, i), X2(:, i));
%
%   dists = cityblkdist(X1, X2, w);
%       computes with weighted city-block distances defined as
%
%           d(x, y) = sum_i w(i) * | x(i) - y(i) |
%
%       w should be an d x 1 column vector specifying the weights of
%       different components, all weights should be non-negative.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Jul 22, 2010
%           - simplify the error handling
%

%% parse and verify input arguments

if ~(ndims(X1) == 2 && isreal(X1) && ndims(X2) == 2 && isreal(X2))
    error('cityblkdist:invalidarg', ...
        'X1 and X2 should be real matrices.');
end

if nargin < 3
    weighted = false;
else
    if ~(isreal(w) && isvector(w))
        error('cityblkdist:invalidarg', ...
            'w should be a real vector.');
    end
    
    if size(w, 1) > 1
        w = w.';
    end
    
    weighted = true;
end


%% main


if ~weighted   
    dists = sum(abs(X1 - X2), 1);
else
    dists = w * abs(X1 - X2);
end

