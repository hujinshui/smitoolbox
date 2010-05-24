function [G, U] = uniqgroup(A, op)
%UNIQGROUP Groups the elements in A based on their values
%
%   [G, U] = uniqgroup(A);
%       groups the elements in A, such that each group of elements take
%       the same value, while different groups have different values.
%
%       The input array A can be anything type admitted by the built-in
%       unique function, including numeric array or cell array of strings.
%
%       In the output, U is the array of distinct values, which equals
%       unique(A); while G is a cell array, with each element being
%       a vector of indices, A(G{i}) is the array of elements whose 
%       values all equal to the i-th element of U.
%
%   [G, U] = uniqgroup(A, 'rows');
%       groups the rows in A, such that each group of rows are the same.
%       U gives the unique rows of A, G gives the corresponding row
%       indices.
%

%   History
%       - Created by Dahua Lin, on May 30, 2008
%

%% main

if nargin == 1
    [U, ~, J] = unique(A);
elseif isequal(op, 'rows')
    [U, ~, J] = unique(A, 'rows');
else
    error('uniqgroup:invalidarg', 'Invalid input arguments.');
end

% group J
[J_sorted, si] = sort(J);
[sp, ep] = valueseg(J_sorted);

% output
G = arrayfun(@(s, e) si(s:e), sp, ep, 'UniformOutput', false);

