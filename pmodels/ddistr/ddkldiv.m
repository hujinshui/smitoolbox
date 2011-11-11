function v = ddkldiv(P, Q)
% Compute Kullback Leibler divergence
%
%   v = ddkldiv(P, Q);
%       computes the Kullback Leibler divergence between the discrete
%       distributions represented by the columns of P and Q.
%
%       If P and Q are both vectors, then it returns the K-L divergence
%       between them.
%
%       If P and Q are matrices of size m x n (with m > 1 and n > 1),
%       then v will be a row vector of size 1 x n, where v(i) is
%       the divergence between P(:,i) and Q(:,i).
%
%   Remarks
%   -------
%       - The value of K-L divergence can be infinity when p_i > 0 and
%         q_i = 0 for some i.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 14, 2010
%       - Modified by Dahua Lin, on Mar 28, 2011
%


%% main

v = sum_xlogy(P, P) - sum_xlogy(P, Q);
   

