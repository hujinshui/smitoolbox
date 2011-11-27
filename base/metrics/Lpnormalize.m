function Y = Lpnormalize(X, p, d)
% Normalize vectors w.r.t. Lp-norm
%
%   Y = Lpnormalize(X, p);
%   Y = Lpnormalize(X, p, d);
%       normalize vectors in X along dimension d w.r.t. Lp-norm. When d is
%       omitted, by default, it is set to 1.
%
%   Y = Lpnormalize(X, p, []);
%       normalize X as a whole w.r.t. Lp-norm.
%

% Created by Dahua Lin, on Aug 1, 2010
%

if nargin < 3
    d = 1;
end

Y = bsxfun(@times, X, 1 ./ Lpnorm(X, p, d));