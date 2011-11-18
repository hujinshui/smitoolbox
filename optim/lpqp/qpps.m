% The function to solve a simple QP on probability simplex
%
%   The problem is formulated as follows:
%
%       min (1/2) ||x||^2 - f' * x;
%
%       s.t. x >= 0 (element wise), and sum(x) == 1.
%
%
%   X = qpps(F);
%
%       When F is a vector, it returns the solution vector.
%
%       When F is an m x n matrix (m > 1 and n > 1), X will be a matrix
%       of the same size. In particular, X(:,i) is the solution
%       corresponding to F(:,i).
%

%   History
%   -------
%       - Created by Dahua Lin, on Mar 28, 2011
%

