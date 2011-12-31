function K = sigmoidkernel(X, Y, a, b)
%SIGMOIDKERNEL Sigmoid kernel
%
%   K = SIGMOIDKERNEL(X, [], a, d);
%   K = SIGMOIDKERNEL(X, Y, a, d);
%
%       Evaluates the sigmoid kernel in a pairwise manner. 
%       
%       The polynomial kernel is defined to be
%
%           k(x, y) = tanh(a * (x' * y) + b)
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
    P = X' * X;
else
    P = X' * Y;
end

K = tanh(a * P + b);

