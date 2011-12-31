function K = polykernel(X, Y, a, b, d)
%POLYKERNEL The polynomial kernel
%
%   K = POLYKERNEL(X, [], a, b, d);
%   K = POLYKERNEL(X, Y, a, b, d);
%
%       Evaluates the polynomial kernel in a pairwise manner. 
%       
%       The polynomial kernel is defined to be
%
%           k(x, y) = (a * (x' * y) + b)^d.
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
K = (a * P + b) .^ d;

