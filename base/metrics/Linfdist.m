function D = Linfdist(X1, X2, d)
%Compute the L_infinity norm distances between corresponding vectors
%
%   D = Linfdist(X1, X2);
%   D = Linfdist(X1, X2, d);
%       computes the L_infinity norm distances between corresponding 
%       vectors in X1 and X2 along dimension d. 
%
%       In the input X1 and X2 should be arrays of the same size.
%
%   D = Linfdist(X1, X2, []);
%       computes the L_infinity norm distances between X1 and X2 as a 
%       whole.
%

% Created by Dahua Lin, on Aug 2, 2010
%

%% main

if nargin < 3
    d = 1;
end

D = Linfnorm(X1 - X2, d);
