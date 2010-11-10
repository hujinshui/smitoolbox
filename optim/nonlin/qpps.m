% Solve a simple quadratic programming problem within probability simplex
%
%   x = qpps(f);
%       This function solves the following quadratic programming problem:
%
%           minimize (1/2) * x^T * x - f^T * x
%
%               s.t. 1^T * x = 1, and x >= 0 (elementwise)
%
%       f can be either a column vector or a d x n matrix, then x
%       will be a column vector or a d x n matrix. In particular, 
%       x(:,i) is the solution corresponding to f(:,i).
%

%   History
%   -------
%       - Created by Dahua Lin, on July 21, 2010
%


