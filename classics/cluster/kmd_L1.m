function R = kmd_L1(X, Y, op)
% K-means related computation (for L1 (cityblock) distance)
%
%   R = kmd_L1(X, Y, 'd');
%   R = kmd_L1(X, Y, 'c');
%       compute the L1 distances between columns in X and columns 
%       in Y pairwisely. 
%
%       Suppose X is a matrix of size d x m and Y is of size d x n,
%       then R is a matrix of size m x n, with R(i, j) being the L1
%       distance between X(:,i) and Y(:,j).
%
%   R = kmd_sqL2(X, [], 'm');
%`      computes the (component-wise) median of X. 
%
%   R = kmd_sqL2(D, [], 't');
%       transforms distance to costs. (Here, it just returns the distance)
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
        R = sum(abs(D), 1);
    else        
        R = zeros(m, size(Y,2));
        for k = 1 : m
            D = bsxfun(@minus, Y, X(:,k));
            R(k,:) = sum(abs(D), 1);
        end
    end
        
elseif op == 'm'
        
    if m > 1
        R = median(X, 2);
    else
        R = X;
    end
    
elseif op == 't'
        
    R = X;    
end


