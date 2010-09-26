function v = L2norm(X, d)
%Compute squared L2-norm values of vectors
%
%   v = L2norm(X);
%   v = L2norm(X, d);
%       compute the L2-norm values of the vectors contained in X along
%       dimension d. If d is omitted, by default, it is set to 1.
%
%   v = L2norm(X, []);
%       compute the L2-norm of the entire array X.
%

% Created by Dahua Lin, on Aug 1, 2010
%

if nargin < 2
    d = 1;
end

v = sqrt(sqL2norm(X, d));