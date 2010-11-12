function [ha, Ja] = gcondupdate(A, X, isigma, w)
% Compute the posterior update from Gaussian conditional distribution
%
%   This function is based on the conditional distribution formulated
%   as follows:
%
%      x ~ N(u, sigma);     (1)
%      x ~ N(A u, sigma);   (2)
%
%   Here, A is a known linear transform matrix. This function computes
%   the posterior update to the Gaussian prior of u, given observed
%   samples x.
%
%   [ha, Ja] = gcondupdate(A, X, isigma);
%   [ha, Ja] = gcondupdate(A, X, isigma, w);
%       computes the posterior update as introduced above. Here, A
%       is the transform matrix. For formulation (1), A should be input
%       as an empty array.
%
%       The observed samples are input as columns of X. isigma is the 
%       inverse of the covariance. 
%
%       The samples can be weighted, and the weights can be input as
%       w. Here, w can be either a row vector of size 1 x n, or a 
%       scalar. If w is omitted, all samples are with a weight 1.
%
%       In the output, ha and Ja are respectively the updates to 
%       potential vectors and information matrix. The formulas to 
%       compute c1a and c2a are given below:
%
%       ha = sum_i w_i * A' * isigma * x_i;
%       Ja = sum_i w_i * A' * isigma * A;
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 17, 2010
%


%% verify input

if ~(isfloat(X) && ndims(X) == 2)
    error('gcondupdate:invalidarg', 'X should be a numeric matrix.');
end

n = size(X, 2);

if ~isempty(A)
    if ~(isfloat(A) && ndims(A) == 2)
        error('gcondupdate:invalidarg', 'A should be a numeric matrix.');
    end
    if size(A, 1) ~= size(X, 1)
        error('gcondupdate:invalidarg', 'The dimension of A and X are inconsistent.');
    end
end

if ~isobject(isigma)
    error('gcondupdate:invalidarg', 'isigma should be a symmetric matrix object.');
end

if nargin < 4
    w = 1;
else
    if ~(isfloat(w) && ...
            (isscalar(w) || (ndims(w) == 2 && size(w,1) == 1 && size(w,2) == n)))
        error('gcondupdate:invalidarg', ...
            'w should be either a scalar or a row vector of length n.');
    end
end


%% main

if isscalar(w)
    if w == 1
        ha = isigma * sum(X, 2);
        Ja = isigma .* n;
    else
        ha = isigma * (sum(X, 2) * w);
        Ja = isigma .* (n * w);
    end
else
    ha = isigma * (X * w');
    Ja = isigma .* sum(w);
end


if ~isempty(A)    
    ha = A' * ha;
    Ja = gsymat(A' * (Ja * A));
end
   

