function r = Linfdiff(A, B)
% Compute Linf-norm of the difference between two arrays
%
%   r = Linfdiff(A, B);
%       computes the L_infinity norm (maximum of absolute value) of the
%       difference between A and B. 
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 12, 2010
%

%% main

D = A - B;
D = abs(D(:));

r = max(D);

