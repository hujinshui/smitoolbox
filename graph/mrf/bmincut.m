function [L, objv] = bmincut(W, h0, h1)
% Solve binary minimum cut problem
%
%   The binary minimum cut problem is formalized as follows:
%
%       minimize sum_{ij} w_{ij} |L_i - L_j| / 2 
%              + sum_i h0(i) |L_i - 0|
%              + sum_i h1(i) |L_i - 1|
%
%   Here, i = 1, ..., n, and each L_i can take values 0 or 1.
%
%   L = bmincut(W, h0, []);
%   L = bmincut(W, [], h1);
%   L = bmincut(W, h0, h1);
%       solves the binary minimum cut problem formalized above and 
%       returns the labeling vector L.
%
%       The inter-node link weights are input as an n x n matrix W.
%       h0 and h1 are row vectors of length n.
%              
%       Either h0 or h1 can be input as an empty array when the
%       corresponding weights are all zeros.
%
%   [L, objv] = bmincut(...);
%       additionally returns the objective value;
%

%% verify input arguments

if ~(isfloat(W) && ndims(W) == 2 && size(W,1) == size(W,2))
    error('bmincut:invalidarg', 'W should be a numeric square matrix.');
end

n = size(W,1);

if ~isempty(h0)
    if ~(isfloat(h0) && ndims(h0) == 2 && size(h0,1) == 1 && size(h0,2) == n)
        error('bmincut:invalidarg', 'h0 should be a numeric vector of size 1 x n.');
    end
end
if ~isempty(h1)
    if ~(isfloat(h1) && ndims(h1) == 2 && size(h1,1) == 1 && size(h1,2) == n)
        error('bmincut:invalidarg', 'h1 should be a numeric vector of size 1 x n.');
    end
end

%% main

% compute h  (h0 - h1)

if ~isempty(h0)
    if ~isempty(h1)
        h = h0 - h1;
    else
        h = h0;
        h1 = 0;
    end
else
    if ~isempty(h1)
        h = -h1;
        h0 = 0;
    else
        error('bmincut:invalidarg', 'Both h0 and h1 are empty.');
    end
end
        
% get edges

[i, j, nw] = find(W);

edges = [i j].';
nw = (nw * 0.5).';

% solve

[L, mf] = graphcut(n, edges, h, nw);

% compute objv

if nargout >= 2
    if isempty(h1)
        objv = mf;
    else
        objv = mf + sum(min(h0, h1));
    end
end

    




