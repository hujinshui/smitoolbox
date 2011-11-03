function D = pwhamdist(X1, X2)
% Compute pairwise Hamming distances
%
%   D = pwhamdist(X1, X2);
%       compute pairwise Hamming distances between the column vectors in
%       X1 and X2. The hamming distance is defined to be the number of
%       different entries in two vectors.
%
%       Let X1 and X2 be d x m and d x n numerical or logical matrices.
%       Then in the output, D will be a matrix of size m x n, where
%       D(i, j) is the hamming distance between X1(:,i) and X2(:,j).
%       
%       Note that this function supports the following value types:
%       double, single, logical, int32, and uint32.
%

%   Created by Dahua Lin, on Aug 2, 2010
%

%% verify input

if nargin < 2 || isempty(X2)
    X2 = X1;
end
    
if ~((isnumeric(X1) || islogical(X1)) && ndims(X1) == 2 && ...
     (isnumeric(X2) || islogical(X2)) && ndims(X2) == 2)
    error('pwhamdist:invalidarg', ...
        'X1 and X2 should be both numeric matrices');
end

if size(X1,1) ~= size(X2,1)
    error('pwhamdist:invalidarg', ...
        'X1 and X2 should have the same number of rows.');
end

%% main

n1 = size(X1, 2);
n2 = size(X2, 2);

D = zeros(n1, n2);

if n1 < n2
    for i = 1 : n1
        D(i, :) = sum(bsxfun(@ne, X2, X1(:,i)), 1) ;
    end
else
    for i = 1 : n2
        D(:, i) = sum(bsxfun(@ne, X1, X2(:,i)), 1).';
    end
end


