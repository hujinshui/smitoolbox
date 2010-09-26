function r = sqL2diff(A, B)
% Compute sqaured L2-norm of the difference between two arrays
%
%   r = L2diff(A, B);
%       computes the squared L2-norm (sum of squared value) of the
%       difference between A and B. 
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 12, 2010
%

%% main

D = A - B;
D = D(:);

r = sum(D.^2);

