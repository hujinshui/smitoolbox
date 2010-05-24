function y = pieces(x, edges, dir, op2)
% PIECES Locates the piece that the input values are in
%
%   y = pieces(x, edges);
%   y = pieces(x, edges, 'left');
%       returns the index of the pieces that each value in x is in.
%       The pieces are delimited by the values in edges.
%       x is assumed to be in non-decreasing order.
%       
%       Suppose, you have n consecutive pieces given by
%       [e_0, e_1), [e_1, e_2), [e_2, e_3), ... [e_{n-1}, e_n),
%       Then the edges should be a vector of length n+1 as
%       [e_0, e_1, ..., e_n].
%
%       In the output, y is of the same size as x, with y(i) giving
%       the piece index for x(i).
%       Concretely, if x(i) has e_{k-1} <= x(i) < e_k, then y(i) = k.
%       If x(i) < e_0, then y(i) = 0, and if x(i) >= e_n, then y(i) = n+1.
%
%   y = pieces(x, edges, 'right');
%       Under this syntax, the pieces are left-open and right-close, as 
%       (e_0, e_1], (e_1, e_2], ..., (e_{n-1}, e_n].
%
%   y = pieces(x, edges, direction, 'unsorted');
%       By default, x is assumed to be in non-decreasing order.
%       If x is not in such order, then one can use the option 'unsorted'.
%
%       In this case, the function will sort the values first, and 
%       re-arrange the output to the original order finally. This 
%       would incur overhead of time-complexity O(nlog(n)).
%   

% Created by Dahua Lin, on Mar 22, 2010
%

%% parse and verify input arguments

assert(isfloat(x) && isvector(x), 'pieces:invalidarg', ...
    'x should be a numeric vector.');

assert(isvector(edges) && isa(edges, class(x)), 'pieces:invalidarg', ...
    'edges should be a numeric vector of the same type as x.');

if nargin < 3 || isempty(dir)
    is_left = true;
else
    assert(ischar(dir), 'pieces:invalidarg', ...
        'The 3rd argument can only be either ''left'' or ''right''.');
    
    if strcmp(dir, 'left')
        is_left = true;
    elseif strcmp(dir, 'right')
        is_left = false;
    else
        error('pieces:invalidarg', ...
            'The 3rd argument can only be either ''left'' or ''right''.');
    end
end        

if nargin < 4
    is_sorted = true;
else
    assert(ischar(op2) && strcmp(op2, 'unsorted'), ...
        'pieces:invalidarg', ...
        'The 4th argument can only be ''unsorted'' if specified.');
    is_sorted = false;
end

%% main

if ~is_sorted
    [x, si] = sort(x);
end

y = pieces_cimp(x, edges, is_left);
y = double(y);

if size(x, 1) > 1
    y = y.';
end

if ~is_sorted
    y(si) = y;
end


