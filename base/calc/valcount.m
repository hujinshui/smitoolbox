function [U, C] = valcount(V)
%Count the number of each unique value in an array
%
%   [U, C] = valcount(V);
%       counts the number of each unique value in the numeric array V.
%
%       In the output, U is a row vector of unique values, C is a 
%       row vector of the same size, where C(i) is the number of 
%       occurrences of U(i) in V.
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

sV = sort(V);
si = valueseg(sV);

U = sV(si);
C = [si(2:end), n+1] - si;

