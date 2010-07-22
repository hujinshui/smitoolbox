function dists = sqeucdist(X1, X2, w)
% Compute squared Euclidean distances between corresponding vectors
%
%   dists = sqeucdist(X1, X2);
%       compute squared Euclidean distances between the corresponding
%       columns of X1 and X2.
%
%       Suppose both X1 and X2 are matrices of size d x n, then dists
%       will be a 1 x n row vector, with dists(i) being the squared
%       Euclidean distance between X1(:,i) and X2(:,i).
%
%   dists = sqeucdist(X1, X2, w);
%       compute weighted squared Euclidean distances between the
%       corresponding vectors of X1 and X2.
%
%       In particular, dists(i) is computed as
%
%           \sum_{j=1}^d w(j) * (x1(j) - x2(j))^2
%
%       Here, x1 and x2 are X1(:,i) and X2(:,i).
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% verify input arguments

if ~(ndims(X1) == 2 && isreal(X1) && ndims(X2) == 2 && isreal(X2))
    error('sqeucdist:invalidarg', ...
        'X1 and X2 should be real matrices.');
end

if nargin < 3
    weighted = false;
else
    if ~(isreal(w) && isvector(w))
        error('sqeucdist:invalidarg', ...
            'w should be a real vector.');
    end
    
    if size(w, 1) > 1
        w = w.';
    end
    
    weighted = true;
end
    


%% compute

D = X1 - X2;

if ~weighted
    dists = sum(D.^2, 1);
else
    D = X1 - X2;    
    dists = w * (D.^2);
end


