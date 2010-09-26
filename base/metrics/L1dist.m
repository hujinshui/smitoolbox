function D = L1dist(X1, X2, d)
%Compute the L1-norm distances between corresponding vectors
%
%   D = L1dist(X1, X2);
%   D = L1dist(X1, X2, d);
%       computes the L1-norm distances between corresponding vectors in 
%       X1 and X2 along dimension d. 
%
%       In the input X1 and X2 should be arrays of the same size.
%
%   D = L1dist(X1, X2, []);
%       computes the L1-norm distances between X1 and X2 as a whole.
%

% Created by Dahua Lin, on Aug 2, 2010
%

%% main

if nargin < 3
    d = 1;
end

D = L1norm(X1 - X2, d);


