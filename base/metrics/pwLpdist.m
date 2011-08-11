function D = pwLpdist(X1, X2, p)
% Compute the pairwise Lp-norm distances
%
%   D = pwLpdist(X1, X2, p);
%       computes the Lp-norm distance between pairs of column vectors in 
%       X1 and X2. 
%
%       Suppose the vector dimension is d, then X1 and X2 should be
%       matrices of size d x m and d x n. In this case, the output
%       is a matrix of size m x n, where D(i, j) is the distance between
%       X1(:,i) and X2(:,j).
%

%   - Created by Dahua Lin, on Aug 2, 2010
%   - Modified by Dahua Lin, on Aug 10, 2011
%       - use pure matlab implementation
%

%% verify input

if nargin < 2 || isempty(X2)
    X2 = X1;
end

if ~(isfloat(X1) && ndims(X1) == 2 && isfloat(X2) && ndims(X2) == 2)
    error('pwLpdist:invalidarg', 'X1 and X2 should be both real matrices.');
end

if size(X1,1) ~= size(X2,1)
    error('pwLpdist:invalidarg', 'X1 and X2 should have the same number of rows.');
end


%% main

if p == 1
    D = pwL2dist(X1, X2);
    
elseif p == 2
    D = pwL2dist(X1, X2);
    
elseif p == inf
    D = pwLinfdist(X1, X2);
    
else
    n1 = size(X1, 2);
    n2 = size(X2, 2);
    
    D = zeros(n1, n2);
    
    if n1 < n2
        for i = 1 : n1
            D(i, :) = sum(abs(bsxfun(@minus, X2, X1(:,i))) .^ p, 1) .^ (1/p);
        end
    else
        for i = 1 : n2
            D(:, i) = sum(abs(bsxfun(@minus, X1, X2(:,i))) .^ p, 1).' .^ (1/p);
        end
    end
end


