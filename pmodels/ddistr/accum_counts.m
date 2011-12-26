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

%% verify input

if isnumeric(V)
    if ~(isfloat(V) && isreal(V) && ndims(V) == 2)
        error('accum_counts:invalidarg', 'V should be a real matrix.');
    end
    uhist = 1;
    [K, n] = size(V);
    
elseif iscell(V) && numel(V) == 2
    K = V{1};
    z = V{2};
    
    if ~(isscalar(K) && isreal(K) && K >= 1 && K == fix(K))
        error('accum_counts:invalidarg', 'K should be a positive integer.');
    end
    
    if ~(isnumeric(z) && isreal(z) && ndims(z) == 2 && size(z, 1) == 1)
        error('accum_counts:invalidarg', 'z should be a numeric row vector.');
    end
    
    uhist = 0;
    n = size(z, 2);    
end

if nargin >= 2 && ~isempty(w)
    if ~(isfloat(w) && isreal(w) && ndims(w) == 2 && size(w, 2) == n)
        error('accum_counts:invalidarg', ...
            'w should be a real matrix with size(w, 2) == n.');
    end
    has_w = 1;
else
    has_w = 0;
end


%% main

if uhist
    
    if has_w
        H = V * w';
    else
        if n == 1
            H = V;
        else
            H = sum(V, 2);
        end
    end
    
else
    
    if has_w
        H = aggreg(w, K, z, 'sum');
    else
        H = intcount(K, z).';
    end
    
end


