function dists = hamdist(X1, X2, w)
%HAMDIST Computes the hamming distances between corresponding vectors.
%
%   dists = hamdist(X1, X2)
%       computes the hamming distances between corresponding vectors in X1
%       and X2.
%
%       Hamming distance between two sequence of 0-1 codes is defined as
%
%           d(x, y) = sum_i 1(x(i) <> y(i))
%
%       which is the number of different code symbols in corresponding
%       positions.
%
%       X1 and X2 should be logical or numeric matrices with the same 
%       size. Let the size be d x n, then dists will be a 1 x n row 
%       vector, with dists(i) being the hamming distance between X1(:,i) 
%       and X2(:,i).
%
%   dists = hamdist(X1, X2, w);
%       computes the weighted hamming distances between corresponding
%       vectors in X1 and X2. 
%
%       The weighted hamming distance is defined as
%
%           d(x, y) = sum_i w_i * 1(x(i) <> y(i))
%
%       The size of w should be d x m, where each column of w gives a 
%       set of weights. The output dists will be a matrix of size 
%       m x n.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Jul 22, 2010
%           - support weighted hamming distances
%           - simplify error handling
%       - Modified by Dahua Lin, on Aug 2, 2010
%           - allow multiple sets of weights
%

%% verify input arguments

if ~(ndims(X1) == 2 && ndims(X2) == 2)
    error('hamdist:invalidarg', ...
        'X1 and X2 should be both matrices.');
end

weighted = (nargin >= 3);


%% main

if ~weighted
    dists = sum(X1 ~= X2, 1);
else
    dists = w.' * (X1 ~= X2);
end

