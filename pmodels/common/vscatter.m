function R = vscatter(X, U, W, op)
% Compute the scatter vector/matrix
%
%   R = vscatter(X, U, W, 'v');
%       computes the scatter vector defined as follows.
%
%       Suppose X is a d x n matrix, U is a d x K matrix, then 
%       R is a d x K matrix defined by
%
%       R(i, k) = sum_j W(k, j) * (X(i, j) - U(i, k))^2.
%
%       Here, W is a K x n matrix. If all weights are equal, the input
%       W can be a scalar.
%
%   R = vscatter(X, U, W, 'c');
%       computes the scatter matrix defined as follows.
%
%       Suppose X is a d x n matrix, U is a d x K matrix, then
%       R is a d x d x K array, defined by
%
%       R(i,i',k) = sum_j W(k,j) * (X(i,j) - U(i,k)) * (X(i',j) - U(i',k))
%
%       Here, W is a K x n matrix. If all weights are equal, the input
%       W can be a scalar.
%
%   Note, U can be input as simply a zero (0), when all its entries are 0.
%       

% Created by Dahua Lin, on Sep 28, 2011
%

%% verify inputs

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('vscatter:invalidarg', 'X should be a real matrix.');
end

if ~(isfloat(U) && isreal(U) && ndims(U) == 2)
    error('vscatter:invalidarg', 'U should be a real matrix.');
end

[d, n] = size(X);

if isequal(U, 0)
    z_u = true;
    K = 1;
else
    z_u = false;
    [d2, K] = size(U);
    if d ~= d2
        error('vscatter:invalidarg', 'The dimensions of X and U are inconsistent.');
    end
end

if ~(isfloat(W) && isreal(W) && (isscalar(W) || isequal(size(W), [K n])))
    error('vscatter:invalidarg', ...
        'W should be either a scalar or a matrix of size K x n.');
end

if ~(ischar(op) && isscalar(op) && (op == 'v' || op == 'c'))
    error('vscatter:invalidarg', ...
        'The 4th argument should be either ''v'' or ''c''.');
end

%% main

if K == 1
    if z_u
        Y = X;
    else
        Y = bsxfun(@minus, X, U);
    end    
    
    if op == 'v'
        if isscalar(W)
            R = sum(Y .^ 2, 2);
            if W ~= 1
                R = R * W;
            end
        else
            R = (Y .^ 2) * W'; 
        end
        
    else % op == 'c'
        if isscalar(W)
            R = Y * Y';
            if W ~= 1
                R = R * W;
            end
        else
            R = Y * bsxfun(@times, Y, W)';
            R = 0.5 * (R + R');
        end
    end     
    
else % K > 1
    
    if op == 'v'
        if isscalar(W)
            T1 = sum(X.^2, 2);
            T2 = bsxfun(@times, sum(X,2), U);
            T3 = n * (U.^2);
            
            R = bsxfun(@plus, T1, T3 - 2 * T2);
            if W ~= 1
                R = R * W;
            end
        else
            T1 = (X.^2) * W';
            T2 = (X * W') .* U;
            T3 = bsxfun(@times, sum(W, 2)', U.^2);
            R = T1 - 2 * T2 + T3;            
        end
        
    else % op == 'c'        
            
        R = zeros(d, d, K);
        
        if isscalar(W)
            for k = 1 : K
                Y = bsxfun(@minus, X, U(:,k));
                C = (Y * Y') * W;                
                R(:,:,k) = C;
            end
        else        
            for k = 1 : K
                Y = bsxfun(@minus, X, U(:,k));
                C = Y * bsxfun(@times, Y, W(k,:))';
                R(:,:,k) = 0.5 * (C + C');            
            end
        end
    end
    
end


