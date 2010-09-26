function v = ddentropy(P)
% Compute the entropy of discrete distribution
%
%   v = ddentropy(P);
%       computes the entropies of the discrete distributions represented
%       by columns of P.
%
%       If P is a vector (row or column), it returns the entropy of 
%       the discrete distribution whose probability mass function is
%       represented by P.
%
%       If P is an m x n (non-vector) matrix,  then v will be a 1 x n 
%       row vector, with v(i) corresponding to P(:,i).
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 18, 2010
%

%% verify 

if ~(isfloat(P) && ndims(P) == 2 && ~issparse(P))
    error('ddentropy:invalidarg', 'P should be a non-sparse numeric matrix.');
end

%% main

v = - safedot(P, log(P));

