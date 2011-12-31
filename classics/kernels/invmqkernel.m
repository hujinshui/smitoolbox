function invmqkernel(X, Y, a, b)
%INVMQKERNEL The inverse multiquadric kernel
%
%   K = INVMQKERNEL(X, [], a, b);
%   K = INVMQKERNEL(X, Y, a, b);
%
%       Evaluates the inverse multi-quadric kernel in a pairwise manner. 
%       
%       The inverse multi-quadric kernel is defined to be
%
%           k(x, y) = 1 / (a * ||x - y||^2 + b)
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

K = 1 ./ (a * D + b);

