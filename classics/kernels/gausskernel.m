function K = gausskernel(X, Y, sigma)
%GAUSSKERNEL The Gaussian kernel (radial basis function)
%
%   K = GAUSSKERNEL(X, [], sigma);
%   K = GAUSSKERNEL(X, Y, sigma);
%
%       Evaluates the Gaussian kernel in a pairwise manner. 
%       
%       The Gaussian kernel is defined to be
%
%           k(x, y) = exp( - ||x - y||^2 / (2 * sigma^2) ).
%
%       Here, each column of X and Y is a sample. Suppose X has m columns
%       and Y has n columns, then K is a matrix of size m x n.
%
%       If Y is input as empty, it means that X and Y are the same.
%       

% Created by Dahua Lin, on Dec 31, 2011
%

%% main

if isempty(Y)
    sx = sum(X.^2, 1);
    D = bsxfun(@minus, bsxfun(@plus, sx.', sx), 2 * (X' * X));
else
    sx = sum(X.^2, 1);
    sy = sum(Y.^2, 1);
    D = bsxfun(@minus, bsxfun(@plus, sx.', sy), 2 * (X' * Y));
end

K = exp(- D * (1 / (2 * sigma^2)));    
    
