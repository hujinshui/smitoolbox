function A = constmat(m, n, v)
% Create a matrix whose entries are fixed to a constant value
%
%   A = constmat(m, n, v);
%       creates a matrix of size m x n, whose entries are fixed to 
%       a constant value v.
%
%       The class of matrix A is is the same as the class of v,
%       which can be either numeric or logical.
%
%   Example:
%
%       A = constmat(2, 3, 10);
%           
%       it returns a matrix A as [10 10 10; 10 10 10]
%

% History
% -------
%   - Created by Dahua Lin, on Jun 5, 2010
%

%% verify input

if ~isscalar(v)
    error('constmat:invalidarg', 'v should be a scalar value.');
end

%% main

if isnumeric(v)
    A = zeros(m, n, class(v));
    A(:) = v;
    
elseif islogical(v)
    if v
        A = true(m, n);
    else
        A = false(m, n);
    end
    
else
    error('constmat:invalidarg', ...
        'v should be either a numeric or a logical value.');
end


