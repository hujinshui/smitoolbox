function [v, g] = logireg_loss(z, y)
%LOGIREG_LOSS The logistic regression loss function
%
%   v = logireg_loss(z, y);
%   [v, g] = logireg_loss(z, y);
%   

% Created by Dahua Lin, on Jan 1, 2012


%% main

pp = 1 ./ (1 + exp(-z));

s0 = y == 0;
s1 = y == 1;

if all(s0 | s1)
    v = zeros(size(z));
    v(s1) = - log(pp(s1));
    v(s0) = - log(1 - pp(s0));
else
    log_pp = log(pp);
    log_pn = log(1 - pp);
    v = - (y .* log_pp + (1-y) .* log_pn);
    v(s1) = - log_pp(s1);
    v(s0) = - log_pn(s0);
end

if nargout >= 2
    g = pp - y;
end

