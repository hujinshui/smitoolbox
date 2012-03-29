function [v, G] = bisquareloss(Z, Y, r)
%BISQUARELOSS Bisquare loss function
%
%   v = BISQUARELOSS(E, r);
%   v = BISQUARELOSS(Z, Y, r);
%
%   [v, G] = BISQUARELOSS(E, r);
%   [v, G] = BISQUARELOSS(Z, Y, r);
%
%       Evaluates the bisquare loss, defined to be
%
%           v = (r^2/6) * min{(1 - (1 - (e/r)^2)^3), 1} 
%
%       If e is a vector, then v is the sum of the loss values at
%       all components. The derivative is given by
%
%           g = e * ( 1 - (e/r)^2 )^2,      when |e| < r
%             = 0                           when |e| >= r
%

% Created by Dahua Lin, on Jan 15, 2012
%

%% main

if nargin == 2
    E = Z;
    r = Y;
elseif nargin == 3
    E = Z - Y;
end

if r == 1
    B = max(1 - E.^2, 0);
else
    B = max(1 - (E * (1/r)).^2, 0);
end

B2 = B.^2;

v = (r^2/6) * (1 - B2 .* B);

if size(v, 1) > 1
    v = sum(v, 1);
end

if nargout >= 2
    G = E .* B2;
end


