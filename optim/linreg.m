function x = linreg(A, y, imat, r)
% Ridge linear regression
%
%   x = linreg(A, y);
%       performs least square linear regression: minimize (1/2) * ||A x- y||^2
%
%       Here, A and y should be matrices of sizes m x n and m x k respectively,
%       then in output x is a matrix of size n x k. 
%
%   x = linreg(A, y, imat);   
%       performs linear regression with a different data norm, i.e.
%       minimize the (1/2) * (Ax - y)' * imat * (Ax - y).
%
%       Here, imat should be an m x m semi-definite matrix, or an
%       m x 1 vector representing a diagonal matrix, or a scalar.
%       Omitting infomat is equivalent to setting it to be 1.
%   
%   x = linreg(A, y, infomat, r);
%       performs ridge regression with the the Tikhonov regularization term
%       given by (1/2) * x' * r * x.
%
%       r can be given in form of either an n x n matrix, an n x 1 vector,
%       or a scalar.
%

%  Created by Dahua Lin, on Nov 26, 2009
%

%% create qp problem

if nargin == 2
    [H, f] = lr2qp(A, y);
    
elseif nargin == 3
    [H, f] = lr2qp(A, y, imat);
    
else
    [H, f] = lr2qp(A, y, imat, r);
    
end

%% solve 

x = - (H \ f);


