function r = L2diff(A, B)
% Compute L2-norm of the difference between two arrays
%
%   r = L2diff(A, B);
%       computes the L2-norm (square root of sum of squared value) of the
%       difference between A and B. 
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 12, 2010
%

%% main

D = A - B;
D = D(:);

r = sqrt(sum(D.^2));

