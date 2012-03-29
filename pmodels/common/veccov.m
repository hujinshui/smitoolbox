function [C, mv] = veccov(X, w, mv)
%VECCOV Computes the covariance matrix of vectors
%
%   C = veccov(X);
%       computes the covariance matrix of column vectors in X.
%
%       The function computes the sample covariance matrix for sample
%       vectors v1, v2, ..., vn, as
%
%           C = sum_i (vi - mv) * (vi - mv)' / n
%
%       where, mv is the mean vector.
%
%       Let X be a d x n matrix with each column representing a sample.
%       Then C will the covariance matrix, whose size is d x d.
%
%   C = veccov(X, w);
%       compute the covariance matrix based on weighted samples. 
%
%       The weighted covariance matrix is computed as
%
%           C = sum_i w(i) * (vi - mv) * (vi - mv)' / (sum_i w(i))
%
%       where, mv is the weighted mean vector.
%
%       Let X be a d x n matrix, then w can be an n x 1 col-vector, with
%       w(i) gives the weight of the sample X(:, i).
%
%       w can also be a n x k matrix, with each col in w giving a set of
%       weights. Then k covariance matrices will be computed, ans thus
%       in the output, C is a d x d x k array, with C(:,:,i) being the
%       covariance matrix computed based on the weights in w(i, :).
%
%   C = veccov(X, [], mv);
%       compute the covariance matrix, with a pre-computed mean vector
%       given by mv. Typically, mv should be a d x 1 vector in a 
%       d-dimensional space.
%
%   C = veccov(X, w, mv);
%       compute the covariance matrix(matrices) using pre-computed 
%       mean vector(s).
%
%       When w is an n x k matrix, mv should be given as a d x k matrix,
%       with mv(:, i) corresponding to the weighted mean vector based on
%       the weights in w(:, i).
%
%   [C, mv] = veccov( ... );
%       additionally returns the mean vector as the 2nd ouptut argument.
%
%   Remarks
%       - Rather than using the original definition given above, the
%         function implements the computation based on the following
%         identity:
%               Cov(X) = E(X * X') - E(X) * E(X)';
%
%         This implementation is more efficient, especially when d < n.
%

%   History
%       - Created by Dahua Lin, on Jun 4, 2008
%       - Modified by Dahua Lin, on Mar 20, 2010
%           - returns mean vector as the second output argument
%       - Modified by Dahua Lin, on April 13, 2010
%

%% parse and verify input arguments

if ~(isfloat(X) && ndims(X) == 2)
    error('veccov:invalidarg', ...
        'X should be a numeric matrix with floating-point value type.');
end

[d, n] = size(X);

if nargin < 2 || isempty(w)    
    weighted = false;
    k = 1;
else    
    if ~(isnumeric(w) && ndims(w) == 2 && size(w,1) == n)
        error('veccov:invalidarg', ...
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
            error('veccov:invalidarg', ...
                'mv should be a d x 1 vector.');
        end
    else
        if ~(isnumeric(mv) && ndims(mv) == 2 && ...
            size(mv,1) == d && size(mv,2) == k)
            error('veccov:invalidarg', ...
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


% compute mean vector(s)
if isempty(mv)
    if ~weighted
        mv = sum(X, 2) / n; 
    else
        mv = X * w;
    end
end

% compute covariance
if ~weighted
    
    C = X * X' * (1/n);    
    if ~all(mv == 0)
        C = C - mv * mv';
    end    
    C = 0.5 * (C + C');
    
else
    if k == 1
        C = X * bsxfun(@times, X', w);
        if ~all(mv == 0)
            C = C - mv * mv';
        end        
        C = 0.5 * (C + C');
    
    else
        C = zeros(d, d, k, class(X));
        
        for i = 1 : k
            cw = w(:, i);
            cmv = mv(:, i);
            
            cc = X * bsxfun(@times, X', cw) - cmv * cmv';
            
            C(:,:,i) = 0.5 * (cc + cc');                
        end
        
    end            
end


