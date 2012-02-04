function [V, mv] = vecvar(X, w, mv)
%VECVAR Variances of vector components
%
%   V = VECVAR(X);
%       computes the variances of the components of vectors in X.
%
%       The variance of values v1, v2, ..., vn is defined as
%
%           (sum_i (v_i - mv)^2) / n
%
%       here, mv is the mean value.
%
%       Given a matrix X of size d x n, with each column representing a
%       d-dimensional vector, then V will be a d x 1 vector, in which
%       V(i) is the variance of the values in the i-th row of X.
%
%   V = VECVAR(X, w);
%       computes the weighted variances of the components of vectors in X.
%
%       The weighted variance is computed as follows
%
%           (sum_i w(i) * (vi - mv)^2) / (sum_i w(i))
%
%       here w(i) is the weight for the i-th value.
%
%       Given a matrix X of size d x n, the weights can be input by a 
%       vector w of size n x 1, with w(i) being the weight for the i-th
%       sample.
%
%       w can also be an n x k matrix, offering multiple sets of different
%       weights. Then, in the output, V will be a d x k matrix, with V(:,i)
%       giving the variances based on the weights in w(i, :).
%
%   V = VECVAR(X, [], mv);
%       computes the weighted variance with pre-computed mean vector. 
%
%       The pre-computed mean vector is given by mv, which should be a
%       d x 1 column vector.
%
%   V = VECVAR(X, w, mv);
%       computes the weighted variances with pre-computed mean vector.
%
%       Then w is a n x 1 column vector, mv should be a d x 1 column vector
%       giving the weighted mean vector based on that weights.
%
%       If w is an n x k matrix offering multiple sets of weights, then
%       mv should be a d x k column vector, with mv(:, i) being the 
%       mean vector based on the weights given in w(:, i).
%
%   [V, mv] = VECVAR( ... );
%       additionally returns the mean vector as the 2nd output argument.
%
%   Remarks
%       - The implementation is based on the following identity rather than
%         the original definition:
%
%           Var(x) = E(x^2) - (E x)^2
%
%         This implementation is typically more efficient, especially when
%         d < n, which is often the case in the context of statistics.
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%       - Modified by Dahua Lin, on Mar 20, 2010
%           - returns mean vector as the second output argument
%       - Modified by Dahua Lin, on April 13, 2010
%

%% parse and verify input arguments

if ~(isfloat(X) && ndims(X) == 2) 
    error('vecvar:invalidarg', ...
        'X should be a numeric matrix with floating-point value type.');
end

[d, n] = size(X);

if nargin < 2 || isempty(w)    
    weighted = false;
    k = 1;
else    
    if ~(isnumeric(w) && ndims(w) == 2 && size(w,1) == n)
        error('vecvar:invalidarg', ...
            'w should be a matrix with n columns.');
    end
    
    k = size(w, 2);
    weighted = true;
end

if nargin < 3
    mv = [];
else
    if k == 1
        if ~(isnumeric(mv) && ndims(mv) == 2 && ...
            size(mv,1) == d && size(mv,2) == 1)
            error('vecvar:invalidarg', ...
                'mv should be a d x 1 vector.');
        end
    else
        if ~(isnumeric(mv) && ndims(mv) == 2 && ...
            size(mv,1) == d && size(mv,2) == k)
            error('vecvar:invalidarg', ...
                'mv should be a d x k matrix.');
        end
    end        
end


%% main

% normalize the weights

if weighted
    if k == 1
        w = w / sum(w);
    else
        w = bsxfun(@times, w, 1 ./ sum(w, 1));
    end
end


% compute E(x)
if isempty(mv)
    if ~weighted
        mv = sum(X, 2) / n;
    else
        mv = X * w;
    end
end

% compute E(x^2)

if ~weighted
    Ex2 = sum(X .* X, 2) / n;
else
    Ex2 = (X .* X) * w;
end

% result

V = Ex2 - mv .* mv;


