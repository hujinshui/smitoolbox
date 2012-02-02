function [v, g] = logistic_loss(z, y)
%LOGISTIC_LOSS The logistic loss function
%
%   v = LOGISTIC_LOSS(z, y);
%   [v, g] = LOGISTIC_LOSS(z, y);
%
%       Evaluates the logistic loss function, which is defined to be
%
%           loss(z, y) = log(1 + exp(-z * y));
%
%       Input arguments:
%       - z:        The linear predictors [1 x n]
%       - y:        The expected response [1 x n]. 
%                   Typically, the values of y can be either 1 or -1.
%       
%       Output arguments:
%       - v:        The loss values [1 x n]
%       - g:        The derivatives [1 x n]
%   

% Created by Dahua Lin, on Jan 1, 2012


%% main

u = 1 + exp(-z .* y);
v = log(u);

if nargout >= 2
    g = (1 - 1 ./ u) .* (-y);
end

