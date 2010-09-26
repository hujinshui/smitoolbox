function Y = L1normalize(X, d)
% Normalize vectors w.r.t. L1-norm
%
%   Y = L1normalize(X);
%   Y = L1normalize(X, d);
%       normalize vectors in X along dimension d w.r.t. L1-norm. When d is
%       omitted, by default, it is set to 1.
%
%   Y = L1normalize(X, []);
%       normalize X as a whole w.r.t. L1-norm.
%

% Created by Dahua Lin, on Aug 1, 2010
%

if nargin < 2
    d = 1;
end

Y = bsxfun(@times, X, 1 ./ L1norm(X, d));