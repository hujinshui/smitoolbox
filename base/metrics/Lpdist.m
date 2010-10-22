function D = Lpdist(X1, X2, p, d)
%Compute the Lp-norm distances between corresponding vectors
%
%   D = Lpdist(X1, X2, p);
%   D = Lpdist(X1, X2, p, d);
%       computes the Lp-norm distances between corresponding vectors in 
%       X1 and X2 along dimension d. 
%
%       In the input X1 and X2 should be arrays of the same size.
%
%   D = Lpdist(X1, X2, p, []);
%       computes the Lp-norm distances between X1 and X2 as a whole.
%

% Created by Dahua Lin, on Aug 2, 2010
%

%% main

if nargin < 4
    d = 1;
end

D = Lpnorm(X1 - X2, p, d);

