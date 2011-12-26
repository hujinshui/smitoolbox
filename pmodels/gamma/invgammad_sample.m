function X = invgammad_sample(alpha, beta, n)
%invgammad_sample Samples from an inverse Gamma distribution
%
%   X = invgammad_sample(alpha, beta);
%   X = invgammad_sample(alpha, beta, n);
%
%       draws n samples from an inverse gamma distribution with shape 
%       parameter alpha and scale parameter beta.
%
%       Input arguments:
%       - alpha:    can be a scalar or a d x 1 column vector.
%       - beta:     can be a scalar or a d x 1 column vector.
%       - n:        the number of samples to draw
%
%       If n is omitted, it is assumed to be 1, and is d is omitted,
%       it is set to size(alpha, 1).
%
%   X = invgammad_sample(alpha, beta, [d, n]);
%       
%       Additionally, specifies the dimension of the samples as d.
%       This syntax is particularly useful when you need to sample
%       from multi-dimensional gamma distributions of which both alpha
%       and beta parameters are scalars.
%

% Created by Dahua Lin, on Sep 1, 2011
% Modified by Dahua Lin, on Dec 26, 2011

%% main

if nargin < 3
    n = 1;
end

X = 1 ./ gammad_sample(alpha, 1 ./ beta, n);

