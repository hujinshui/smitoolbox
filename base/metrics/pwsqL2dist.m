function D = pwsqL2dist(X1, X2, w)
% Compute the pairwise squared L2-norm distances
%
%   D = pwsqL2dist(X1, X2);
%       computes the squared L2-norm distance between pairs of column 
%       vectors in X1 and X2. 
%
%       Suppose the vector dimension is d, then X1 and X2 should be
%       matrices of size d x m and d x n. In this case, the output
%       is a matrix of size m x n, where D(i, j) is the distance between
%       X1(:,i) and X2(:,j).
%
%   D = pwsqL2dist(X);
%   D = pwsqL2dist(X, []);
%       computes the squared L2-norm distance between pairs of column
%       vectors in X. The implementation for this case is more efficient
%       than pwsqL2dist(X, X), despite that both yield the same result.
%
%   D = pwsqL2dist(X, [], w);
%   D = pwsqL2dist(X1, X2, w);
%       computes the weighted squared L2-norm distance between column 
%       vectors in X1 and X2. The weighted squared L2-norm distance is 
%       defined by
%
%           d = sum_i w(i) * |x1(i) - x2(i)|^2
%
%       In the input, w should be a column vector.
%

%   Created by Dahua Lin, on Aug 2, 2010
%

%% verify input

if nargin < 2
    X2 = [];
end

if ndims(X1) ~= 2 || ~isreal(X1)
    error('pwsqL2dist:invalidarg', 'X1 should be a real matrix.');
end

if ~isempty(X2) && (ndims(X2) ~= 2 || ~isreal(X2))
    error('pwsqL2dist:invalidarg', 'X2 should be a matrix.');
end

if nargin < 3 || isempty(w)
    weighted = false;
else
    if ~(ndims(w) == 2 && size(w, 2) == 1 && isreal(w))
        error('pwsqL2dist:invalidarg', 'w should be a real row vector.'); 
    end
    weighted = true;
end


%% main

if isempty(X2)   
        
    n = size(X1, 2);
    
    if ~weighted    
        D = X1' * X1;
    else
        D = X1' * bsxfun(@times, w, X1);
    end
    
    % take the diagonal elements, which equals sum(X1 .* X1, 1);
    sx = D(1 + (n+1) * (0:n-1));
    
    D = bsxfun(@plus, (-2) * D, sx);
    D = bsxfun(@plus, D, sx');        
    
else
    
    if ~weighted    
        D = (-2) * (X1' * X2);            
        D = bsxfun(@plus, D, sum(X1 .^ 2, 1).');
        D = bsxfun(@plus, D, sum(X2 .^ 2, 1));
    else
        wX1 = bsxfun(@times, w, X1);
        wX2 = bsxfun(@times, w, X2);
        D = (-2) * (X1' * wX2);
        D = bsxfun(@plus, D, sum(X1 .* wX1, 1).');
        D = bsxfun(@plus, D, sum(X2 .* wX2, 1));
    end
        
end


% enforce non-negativeness (rounding error may lead to negative values)
D(D < 0) = 0;

