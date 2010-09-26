function v = Linfnorm(X, d)
%Compute L_infinity norm values of vectors
%
%   v = Linfnorm(X);
%   v = Linfnorm(X, d);
%       compute the L_infinity norm values of the vectors contained in X 
%       along dimension d. If d is omitted, by default, it is set to 1.
%
%   v = Linfnorm(X, []);
%       compute the L_infinity norm of the entire array X.
%

% Created by Dahua Lin, on Aug 1, 2010
%

%% verify input

if nargin < 2
    d = 1;
end

%% main

if ~isempty(d)
    v = max(abs(X), [], d);
else
    v = max(abs(X(:)));
end
