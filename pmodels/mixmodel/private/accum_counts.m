function H = accum_counts(V, w)
% Accumulate counts into histograms
%
%   H = accum_counts(V);
%   H = accum_counts(V, w);
%
%       accumulates the (weighted) counts in V.
%
%       Input arguments:
%       - V:    the observed histograms, a matrix of size K x n.
%               Here, K is the number of classes, and n is the number of
%               observed histograms.
%
%               V can also be in form of a cell array as {K, z}, where
%               z is a row vector of size 1 x n, and each element in z
%               is an integer between 1 and K.
%
%       - w:    the weights of the histograms, a row vector of size 1 x n
%               or a matrix of size m x n (m groups of weights).
%
%       Output arguments:
%       - H:    the accumulated histograms
%               if w is omitted or a row vector, then P is a column vector
%               of size K x 1, if there are multiple groups of weights,
%               namely m > 1, then P is of size K x m.
%

% Created by Dahua Lin, on Dec 25, 2011
%


%% main

if nargin < 2
    w = [];
end

if isnumeric(V)
    
    if isempty(w)
        if n == 1
            H = V;
        else
            H = sum(V, 2);
        end
    else
        H = V * w';
    end
    
elseif iscell(V)
    
    K = V{1};
    z = V{2};

    if isempty(w)
        H = intcount(K, z).';
    else
        H = aggreg(w, K, z, 'sum');
    end
    
end


