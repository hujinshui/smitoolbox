% Find the set of indices of the occurrences of integers in a given range
%
%   G = intgroup([v0, v1], V);
%       finds the set of indices of occurrences of integers in range
%       [v0, v1] in the array V.
%
%       In the input, V is a numeric array containing integers (whose
%       value class can be either double, single, or int32).
%       In the output, G is a cell array of size 1 x (v1 - v0 + 1),
%       where G{i} is the indices of entries in V whose values equal
%       v0 + i - 1.
%

%   History
%   -------
%       - Created by Dahua Lin, on June 6, 2010
% 


