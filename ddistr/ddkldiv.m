function v = ddkldiv(P, Q)
% Compute Kullback Leibler divergence
%
%   v = ddkldiv(P, Q);
%       computes the Kullback Leibler divergence between the discrete
%       distributions represented by the columns of P and Q.
%
%       P and Q should be matrices with the same number of rows. Each row 
%       is a vector of non-negative entries that sum to one. The sizes of
%       P and Q can be both m x n, or in such a way that either one is 
%       m x 1, while the other one is m x n. For the latter case, the
%       m x 1 one will be broadcasted. For both cases, the size of v will
%       be 1 x n.
%
%   Remarks
%   -------
%       - The value of K-L divergence can be infinity when p_i > 0 and
%         q_i = 0 for some i.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 14, 2010
%

%% verify input

if ~(isfloat(P) && isfloat(Q) && ndims(P) == 2 && ndims(Q) == 2 && ...
        size(P, 1) == size(Q, 1))
    error('ddkldiv:invalidarg', ...
        'P and Q should be numeric matrices with the same number of rows.');
end

n1 = size(P, 2);
n2 = size(Q, 2);

if ~(n1 == n2 || n1 == 1 || n2 == 1)
    error('ddkldiv:invalidarg', ...
        'The number of columns in P and Q are inconsistent.');
end

%% main

v = safedot(P, log(P)) - safedot(P, log(Q));
   

