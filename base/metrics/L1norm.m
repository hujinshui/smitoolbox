function v = L1norm(X, d)
%Compute L1-norm values of vectors
%
%   v = L1norm(X);
%   v = L1norm(X, d);
%       compute the L1-norm values of the vectors contained in X along
%       dimension d. If d is omitted, by default, it is set to 1.
%
%   v = L1norm(X, []);
%       compute the L1-norm of the entire array X.
%

% Created by Dahua Lin, on Aug 1, 2010
%

%% verify input

if nargin < 2
    d = 1;
end

%% main

if ~isempty(d)
    v = sum(abs(X), d);
else
    v = sum(abs(X(:)));
end

