function [ha, Ja] = gcondupdate(A, X, icov, cf, w)
% Compute the posterior update from Gaussian conditional distribution
%
%   This function is based on the conditional distribution formulated
%   as follows:
%
%      x ~ N(u, C);     (1)
%      x ~ N(A u, C);   (2)
%
%   Here, A is a known linear transform matrix. This function computes
%   the posterior update to the Gaussian prior of u, given observed
%   samples x.
%
%   [ha, Ja] = gcondupdate(A, X, icov, cf);
%   [ha, Ja] = gcondupdate(A, X, icov, cf, w);
%
%       Input arguments:
%       - A:    the transformation matrix as in Eq.(2). When A is an 
%               identity transform, it can be input as an empty array [].
%       - X:    the matrix comprised of observed samples as columns.
%       - icov: the inverse of covariance C.
%       - cf:   the char specifying the form of icov.
%       - w:    the sample weights, which can be either a row vector of
%               size 1 x n, or a scalar.
%               If w is omitted, all samples are with a weight 1.
%
%       In the output, ha and Ja are respectively the updates to 
%       potential vectors and information matrix. The formulas to 
%       compute c1a and c2a are given below:
%
%       When A is [], 
%       ha = sum_i w_i * icov * x_i;
%       Ja = sum_i w_i * icov;
%
%       When A is non-empty, then
%       ha = sum_i w_i * A' * isigma * x_i;
%       Ja = sum_i w_i * A' * isigma * A;
%
%       Note, when A isempty, Ja has the same covariance form as icov,
%       when A is non-empty, Ja is the full covariance matrix.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 17, 2010
%       - Modified by Dahua Lin, on Aug 18, 2011
%


%% verify input

if ~(isempty(A) || (isfloat(A) && isreal(A) && ndims(A) == 2))
    error('gcondupdate:invalidarg', 'A should be either empty or a real matrix.');
end

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('gcondupdate:invalidarg', 'X should be a real matrix.');
end
[d, n] = size(X);

if ~isempty(A)
    if size(A, 1) ~= d
        error('gcondupdate:invalidarg', 'The dimension of A and X are inconsistent.');
    end
end

if ~(ischar(cf) && isscalar(cf))
    error('gcondupdate:invalidarg', 'The covariance form cf should be a char scalar.');
end

switch cf
    case 's'
        icov_ok = (isfloat(icov) && isreal(icov) && isscalar(icov));           
    case 'd'
        icov_ok = (isfloat(icov) && isreal(icov) && isequal(size(icov), [d 1]));        
    case 'f'
        icov_ok = (isfloat(icov) && isreal(icov) && isequal(size(icov), [d d])); 
    otherwise
        error('gcondupdate:invalidarg', 'The covariance form cf is invalid.');
end
if ~icov_ok
    error('gcondupdate:invalidarg', 'The argument icov is invalid.');
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
    ha = gmat_mvmul(cf, icov, sum(X, 2));
    if w ~= 1
        ha = ha * w;
    end
    Ja = icov * (n * w);
else
    ha = gmat_mvmul(cf, icov, (X * w'));
    Ja = icov * sum(w);
end


if ~isempty(A)    
    ha = A' * ha;        
    Ja = A' * gmat_mvmul(Ja * A);
    Ja = 0.5 * (Ja + Ja');
end
   

