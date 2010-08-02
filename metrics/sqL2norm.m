function v = sqL2norm(X, d)
%Compute squared L2-norm values of vectors
%
%   v = sqL2norm(X);
%   v = sqL2norm(X, d);
%       compute the squared L2-norm values of the vectors contained in X 
%       along dimension d. If d is omitted, by default, it is set to 1.
%
%   v = sqL2norm(X, []);
%       compute the squared L2-norm of the entire array X.
%

% Created by Dahua Lin, on Aug 1, 2010
%

%% verify input

if nargin < 2
    d = 1;
end

%% main

if ~isempty(d)
    v = sum(X.^2, d);
else
    v = sum(X(:).^2);
end
