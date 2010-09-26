function Y = Linfnormalize(X, d)
% Normalize vectors w.r.t. L_infinity norm
%
%   Y = Linfnormalize(X);
%   Y = Linfnormalize(X, d);
%       normalize vectors in X along dimension d w.r.t. L_infinity norm. 
%       When d is omitted, by default, it is set to 1.
%
%   Y = Linfnormalize(X, []);
%       normalize X as a whole w.r.t. L_infinity norm.
%

% Created by Dahua Lin, on Aug 1, 2010
%

if nargin < 2
    d = 1;
end

Y = bsxfun(@times, X, 1 ./ Linfnorm(X, d));