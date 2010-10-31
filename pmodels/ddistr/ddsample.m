function S = ddsample(P, n, randStr)
% Draw samples from discrete distribution
%
%   S = ddsample(P);
%   S = ddsample(P, n);
%       If P is a vector, then it draws n samples from the discrete
%       distribution whose probability mass function is given by P.
%       It returns a vector S of length n, which is a row or a column,
%       depending on the shape of P. 
%
%       The values in S are integers in [1, m], where m is the length
%       of P (the number of classes). However, if sum(P) < 1, then
%       it is possible to draw a value m+1, with probability 1 - sum(P).
%
%       If P is a non-vector matrix, it regards each column of P as a
%       discrete distribution, and draws n samples from each.
%       Suppose P has c columns, then S is a matrix of size n x c, 
%       where S(:,i) are the samples corresponding to P(:,i).
%
%       When n is omitted, it is set to 1 by default.
%   
%   S = ddsample(P, n, randStr);
%       Do the sampling with the given random number stream
%

% Created by Dahua Lin, on Sep 13, 2010
%

%% verify input

if ~(ndims(P) == 2 && isfloat(P))
    error('ddsample:invalidarg', 'P should be a numeric matrix.');
end

if nargin < 2
    n = 1;
end

if nargin < 3
    randStr = [];
end
  
%% main

if isvector(P)    
    if size(P, 1) > 1
        F = [0; cumsum(P)];
        if isempty(randStr)
            rv = rand(n, 1);
        else
            rv = randStr.rand(n, 1);
        end
    else
        F = [0, cumsum(P)];
        if isempty(randStr)
            rv = rand(1, n);
        else
            rv = randStr.rand(1, n);
        end
    end
    
    S = pieces(rv, F);        
else
    nc = size(P, 2);
    F = [zeros(1, nc); cumsum(P, 1)];    
    if isempty(randStr)
        rv = rand(n, nc);
    else
        rv = randStr.rand(n, nc);
    end
    
    S = zeros(n, nc);
    for i = 1 : nc
        S(:, i) = pieces(rv(:,i), F(:,i));
    end    
end


