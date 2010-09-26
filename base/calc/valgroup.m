function [U, G] = valgroup(V)
% Find the set of indices for each unique value
%
%   [U, G] = valgroup(V);
%       find the set of indices for each unique value appeared in V.
%
%       For the input, V should be a numeric array. For the output, 
%       U is a row vector of sorted unique values, while G is a cell 
%       array, and G{i} is a row vector of indices whose corresponding
%       values in V are U(i).
%

%   History
%   -------
%       - Created by Dahua Lin, on June 7, 2010
%

%% verify input

if ~isnumeric(V)
    error('valcount:invalidarg', ...
        'V should be a numeric array.');
end

%% main

n = numel(V);
if ndims(V) > 2 || size(V, 2) > 1
    V = reshape(V, 1, numel(V));
end

[sV, si] = sort(V);
sp = valueseg(sV);
ep = [sp(2:end) - 1, n];

U = sV(sp);

m = numel(sp);
G = cell(1, m);
for i = 1 : m
    G{i} = si(sp(i):ep(i));
end

