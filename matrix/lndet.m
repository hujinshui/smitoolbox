function v = lndet(C)
% Computes the logarithm of determinant of a positive definite matrix
%
%   v = lndet(C);
%       returns the log-determinant of the positive definite matrix C.
%
%   Remarks
%   -------
%       This implementation directly computes the log-determinant of
%       C based on Cholesky decomposition, and thus can effectively 
%       avoid the issue of overflow or underflow.
%

L = chol(C);
v = 2 * sum(log(diag(L)));
