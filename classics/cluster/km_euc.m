function R = km_euc(X, Y, op)
% K-means related computation (for Euclidean distance)
%
%   R = km_euc(X, Y, 'd');
%       compute the Euclidean distances between columns in X and columns 
%       in Y pairwisely. 
%
%       Suppose X is a matrix of size d x m and Y is of size d x n,
%       then R is a matrix of size m x n, with R(i, j) being the Euclidean
%       distance between X(:,i) and Y(:,j).
%
%   R = km_euc(X, Y, 'c');
%       compute the matching costs between columns in X and those in Y.
%
%       Suppose X is a matrix of size d x m and Y is of size d x n,
%       then R is a matrix of size m x n, with R(i, j) being the cost of
%       assigning Y(:,j) to the cluster with center X(:,i). Here, it is
%       the squared Euclidean distance.
%
%   R = km_euc(X, [], 'm');
%   R = km_euc(X, w, 'm');
%`      computes the (weighted) mean of X. When w is an empty array, it
%       computes the un-weighted mean.
%
%   R = km_euc(D, [], 't');
%       transforms distance to costs. (Here, it computes squares of D).
%
%   Remarks
%   -------
%       - This function implements the computation for Euclidean distance
%         in K-means. To use other types of distance, one can implement
%         his/her own computation function with the same interface.
%

% Created by Dahua Lin, on Sep 27, 2010
%

%% parse and verify input

m = size(X, 2);

%% main

if op == 'c' || op == 'd'
    
    if m == 1
        D = bsxfun(@minus, Y, X);
        R = sum(D .^2, 1);
    else
        R = bsxfun(@plus, bsxfun(@plus, (-2) * X' * Y, sum(Y.^2, 1)), sum(X.^2, 1)');        
        R = max(R, 0);  % neg values may appear due to finite-precision computation
    end
    
    if op == 'd'
        R = sqrt(R);
    end
    
elseif op == 'm'
    
    w = Y;    
    if m > 1
        if isempty(w)
            R = sum(X, 2) * (1 / size(X, 2));
        else
            R = X * (w' / sum(w));
        end
    else
        R = X;
    end
    
elseif op == 't'
    
    D = X;
    R = D.^2;
end


