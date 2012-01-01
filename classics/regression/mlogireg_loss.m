function [v, g] = mlogireg_loss(Z, Y)
%MLOGIREG_LOSS Multi-class Logistic regression Loss function
%
%   v = mlogireg_loss(Z, Y);
%   [v, g] = mlogireg_loss(Z, Y);
%
%

% Created by Dahua Lin, on Jan 1, 2012
%

%% main

[K, n] = size(Z);
P = nrmexp(Z, 1);

if size(Y, 1) > 1
    v = - sum_xlogy(Y, P, 1);
    
    if nargout >= 2
        g = P - Y;
    end    
else    
    si = Y + (0:n-1) * K; 
    p = P(si);
    v = - log(p);
            
    if nargout >= 2
        g = P;
        g(si) = g(si) - 1;
    end
end    
