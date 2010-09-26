function r = L1diff(A, B)
% Compute L1-norm of the difference between two arrays
%
%   r = L1diff(A, B);
%       computes the L1-norm (sum of absolute value) of the difference
%       between A and B. 
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 12, 2010
%

%% main

D = A - B;
D = abs(D(:));

r = sum(D);

