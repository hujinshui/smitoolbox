function D = L2dist(X1, X2, d)
%Compute the L2-norm distances between corresponding vectors
%
%   D = L2dist(X1, X2);
%   D = L2dist(X1, X2, d);
%       computes the L2-norm distances between corresponding vectors in 
%       X1 and X2 along dimension d. 
%
%       In the input X1 and X2 should be arrays of the same size.
%
%   D = L2dist(X1, X2, []);
%       computes the L2-norm distances between X1 and X2 as a whole.
%

% Created by Dahua Lin, on Aug 2, 2010
%

%% main

if nargin < 3
    d = 1;
end

D = L2norm(X1 - X2, d);

