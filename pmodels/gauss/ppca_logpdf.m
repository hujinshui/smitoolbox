function L = ppca_logpdf(M, X)
% Evaluate log pdf values at given samples
%
%   L = ppca_logpdf(M, X);
%       evaluates the log proability density values at the samples given 
%       by X, with respect to the PPCA model.
%
%       Suppose X has n columns (each column is sample), then L
%       is a vector of size 1 x n.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% main

dists = ppca_sqmahdist(M, X);
L = (-0.5) * ((M.d * log(2*pi) + M.ldc) + dists);
