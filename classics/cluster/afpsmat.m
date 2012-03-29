function S = afpsmat(S0, q)
%AFPSMAT Affinity propagation similarity matrix
%
%   S = AFPSMAT(S0);
%       
%       Given an input similarity matrix, it generates a modified matrix
%       for the use in affinity propagation.
%
%       In particular, it sets the median of the pairwise similarities 
%       to each diagonal entry of S.
%
%   S = AFPSMAT(S0, q);
%
%       It sets the diagonal entries to the q-quantile of the pairwise
%       similarities.
%
%       AFPSMAT(S0, 0.5) is equivalent to AFPSMAT(S0);
%
%
%   Remarks
%   -------
%       if you have a pairwise distance(dissimilarity) matrix D, you
%       can input -D as a similarity matrix.
%

% Created by Dahua Lin, on Mar 28, 2012
%

%% verify input arguments

n = size(S0, 1);
if ~(isfloat(S0) && isreal(S0) && ~issparse(S0) && size(S0, 2) == n)
    error('afpsmat:invalidarg', ...
        'S0 should be a non-sparse real square matrix.');
end

if nargin < 2
    q = 0.5;
else
    if ~(isfloat(q) && isscalar(q) && isreal(q) && q > 0 && q < 1)
        error('afpsmat:invalidarg', ...
            'q should be a real scalar with 0 < q < 1.');
    end
end


%% main

ss = S0;
ss(1:(n+1):n^2) = [];

if q == 0.5
    dv = median(ss);
else
    dv = quantile(ss, q);
end

S = S0;
S(1:(n+1):n^2) = dv;

