function [v, d1, d2] = svm_huber_loss(z, h)
% Compute the Huber loss with margin 1 as follows
%
%   f(z) = 
%       0,              if z >= 1 + h
%       (1+h-z)^2/4h,   if |1 - z| < h
%       1 - z,          if z <= 1 - h
%
%   [v, d1, d2] = svm_huber_loss(z, h);
%
%   Inputs:
%       - z:    the array of input values
%       - h:    the quadratic band-width [scalar]
%               if h = 0, it reduces to hinge loss
%       
%   Outputs:
%       - v:    the function values
%       - d1:   the first derivative values
%       - d2:   the second derivative values
%

% Created by Dahua Lin, on Apr 21, 2011
%

%% main

u = 1 - z;

if h > 0
    sq = find(abs(u) < h);
    s1 = find(u >= h);
else
    s1 = find(u > h);
end

v = zeros(size(u));
if h > 0
    v(sq) = ((u(sq) + h).^2 / (4*h));
end
v(s1) = u(s1);

if nargout >= 1
    d1 = zeros(size(u));
    if h > 0
        d1(sq) = - ((u(sq) + h) / (2*h));
    end
    d1(s1) = -1;
end

if nargout >= 2
    d2 = zeros(size(u));
    if h > 0
        d2(sq) = 1 / (2*h);
    end
end
    
