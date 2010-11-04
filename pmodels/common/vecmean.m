function mv = vecmean(X, w)
%VECMEAN Computes the mean of vectors
%
%   mv = vecmean(X);
%       Computes the mean vector of the column vectors in X.
%
%       The mean vector of vectors v1, v2, ..., vn is defined as
%   
%           mv = (v1 + v2 + ... + vn) / n;
%
%       Let X be a d x n matrix, with each column representing a
%       d-dimensional sample. Then mv is a d x 1 column vector,
%       which is their mean vector.
%
%   mv = vecmean(X, w);
%       Computes the weighted mean vector of the columns in X.
%
%       The weighted mean vector is defined as
%
%           mv = sum_i w(i) * vi / (sum_i w(i))
%
%       Let X be a d x n matrix, then w can be a 1 x n row vector,
%       in which, w(i) is the weight of the i-th sample in X.
%
%       w can also be a k x n matrix, with each row giving a set
%       of weights. In this case, the output mv will be a d x k
%       matrix, with mv(:, i) being the weighted mean vector
%       computed based on the weights in w(i, :).
%

%   History
%       - Created by Dahua Lin, on Jun 4, 2008
%       - Modified by Dahua Lin, on April 13, 2010
%

%% parse and verify input arguments

if ~(isfloat(X) && ndims(X) == 2) 
    error('vecmean:invalidarg', ...
        'X should be a numeric matrix with floating point value type.');
end

n = size(X, 2);

if nargin < 2 || isempty(w)    
    weighted = false;
else    
    if ~(isnumeric(w) && ndims(w) == 2 && size(w,2) == n)
        error('vecmean:invalidarg', ...
            'w should be a matrix with n columns.');
    end    
    weighted = true;
end

%% main

if ~weighted
    mv = sum(X, 2) * (1/n);
else
    if size(w, 1) == 1
        mv = X * (w / sum(w))';
    else
        mv = X * bsxfun(@times, w, 1 ./ sum(w, 2))';
    end
end

    