function n = numuniq(A, op)
%NUMUNIQ Counts the number of unique (distinct) values
%
%   n = numuniq(A);
%       returns the number of distinct values in A.
%
%   n = numuniq(A, 'rows');
%       returns the number of distinct rows in A.
%

%   History
%       - Created by Dahua Lin, on May 31, 2008
%

%% main

if nargin == 1
    n = numel(unique(A));
    
elseif nargin == 2 && isequal(op, 'rows')
    n = size(unique(A), 1);
    
else
    error('numuniq:invalidarg', 'Invalid input arguments to numuniq');
end
    