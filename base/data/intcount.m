function c = intcount(rgn, V, dim)
%INTCOUNT Count the number of occurrence of integers
%
%   c = INTCOUNT([v0, v1], V);
%       counts the number of integers in [v0, v1] that occurs in V.
%       
%       In the input, V is an array containing integer values (whose
%       type can be double or single nonetheless). 
%       In the output, c is a row vector of length v1 - v0 + 1.
%       In particular, c(i) counts the number of times that the value
%       v0 + (i-1) occurs in V.
%
%   c = INTCOUNT(K, V);
%       This is equivalent to intcount([1, K], V);
%
%   c = INTCOUNT([v0, v1], V, dim);
%   c = INTCOUNT(K, V, dim);
%
%       performs the counting along specified direction. Here, dim
%       can be either 1 or 2. 
%
%       If dim == 1, it counts for each column, resulting a vector
%       of counts, whose size is 1 x n.
%
%       If dim == 2, it counts for each row, resulting a vector of 
%       counts, whose size is m x 1.
%

%   History
%   -------
%       - Created by Dahua Lin, on May 26, 2010
%       - Modified by Dahua Lin, on Jun 6, 2010
%           - Use pure C++ mex implementation
%       - Modified by Dahua Lin, on Feb 26, 2012
%           - Support per column/row counting
%

%% verify input arguments

if isnumeric(rgn) && isvector(rgn) 
    if isscalar(rgn)
        v0 = 1;
        v1 = rgn;
        if v1 < 1
            error('intcount:invalidarg', 'K should be a positive number.');
        end        
    elseif numel(rgn) == 2
        v0 = rgn(1);
        v1 = rgn(2);
        if v0 > v1
            error('intcount:invalidarg', 'The condition v0 <= v1 is not met.');
        end        
    else
        error('intcount:invalidarg', 'The first argument is invalid.');
    end
else
    error('intcount:invalidarg', 'The first argument is invalid.');
end
v0 = int32(v0);
v1 = int32(v1);

if ~(isnumeric(V) && isreal(V) && ~issparse(V))
    error('intcount:invalidarg', ...
        'V should be a non-sparse numeric array of integers.');
end

if nargin < 3
    dim = 0;
else
    if ndims(V) > 2
        error('intcount:invalidarg', ...
            'V need to be a matrix when dim is explicitly specified.');
    end    
    if isequal(dim, 1)
        dim = 1;
    elseif isequal(dim, 2)
        dim = 2;
        V = V.';
    end
end

%% main

c = intcount_cimp(v0, v1, V, dim > 0);

if dim == 2
    c = c.';
end

    