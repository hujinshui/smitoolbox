function r = fast_median(X, dim)
% Find the median value using O(n) algorithm
%
%   r = fast_median(X);
%       If X is a vector, it returns the median of X. 
%       If X is a matrix of size m x n where m > 1 and n > 1, then
%       it returns a row vector r, such that r(j) is the median value
%       for X(:, j).
%
%   r = fast_median(X, dim);
%       compute the median values along the dimension specified by dim.
%
%   Remarks
%   -------
%       - Note the current version only supports the case where X is
%         a matrix.
%
%       - The implementation is based on an O(n) selection algorithm
%         in C++ STL (std::nth_element).
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 11, 2010
%       - Modified by Dahua Lin, on Mar 28, 2010
%

%% verify input arguments

if ~(isfloat(X) && ~issparse(X) && isreal(X))
    error('fast_median:invalidarg', 'X should be a non-sparse real array.');
end

if nargin < 2
    dim = 0;
else
    if ~(isnumeric(dim) && isscalar(dim) && (dim == 1 || dim == 2))
        error('fast_median:invalidarg', 'dim should be either 1 or 2');
    end
    dim = double(dim);
end


%% main

if ~isempty(X)    
    r = fast_median_cimp(X, dim);
    
else
    % let median itself to deal with empty cases for conformant behavior
    
    if nargin == 1
        r = median(X);
    else
        r = median(X, dim);
    end    
end





