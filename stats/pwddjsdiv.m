function D = pwddjsdiv(P, Q)
%PWDDJSDIV Computes pairwise J-S divergences
%
%   D = pwjsdiv(P, Q);
%       computes pairwise Jensen-Shannon divergence between the discrete
%       distributions in P and Q.
%
%       Jensen-Shannon divergence, also known as Jeffrey divergence,
%       is a "symmetric version" of the well-known Kullback-Leibler 
%       divergence. J-S divergence between p and q is defined as
%
%       JSD(p || q) = (1/2) sum_i p(i) log p(i)/u(i) + q(i) log q(i)/u(i)
%
%       Here, u = (p + q)/2, which is the "average distribution".
%
%       Compared to K-L divergence, J-S divergence is symmetric and
%       numerically stable.
%
%       Let P and Q be d x m and d x n matrices respectively, then 
%       D will be an m x n matrix, with D(i, j) being the J-S divergence
%       between the probabilities with their probability mass function
%       given in P(:, i) and Q(:, j).
%
%   D = pwddjsdiv(P);
%       computes the pairwise J-S divergence between the discrete
%       distributions in columns of P. It is equivalent to 
%       pwddjsdiv(P, P);           
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'pwddjsdiv:invalidarg', ...
    'P should be a probability mass matrix with all elements within [0, 1]');

if nargin < 2
    Q = [];
else
    assert(isfloat(Q) && ndims(Q) == 2 && all(Q(:) >= 0 & Q(:) <= 1), ...
        'pwddjsdiv:invalidarg', ...
        'Q should be a probability mass matrix with all elements within [0, 1]');
end

assert(size(P,1) == size(Q,1), ...
    'pwddjsdiv:invalidarg', ...
    'P and Q should have the same number of rows.');


%% main

if isempty(Q)
    pe = ddentropy(P);
    qe = pe; 
    
    Q = P;
    
else
    pe = ddentropy(P);
    qe = ddentropy(Q);

end

D = - pwjsdiv_ulogu_cimp(P, Q);

D = bsxfun(@minus, D, 0.5 * pe');
D = bsxfun(@minus, D, 0.5 * qe);

