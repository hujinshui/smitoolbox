function [nbs, nbdists] = findnbs(dists, K, thres)
% Find neighbors of points based on pairwise distances
%
%   nbs = findnbs(dists, K);
%       finds the K-nearest neighbor of each point based on the pairwise
%       distance matrix given in dists.
%
%       Suppose there are n points, dists should be an n x n matrix.
%       In addition, K should be a positive value with K < n.
%
%       In the output, nbs is a cell array, with nbs{i} being a column
%       vector of neighboring sample indices of the i-th element.
%       The indices are sorted in ascending order of distance.
%       
%   nbs = findnbs(dists, K, thres);
%       finds the at most K nearest neighbors of each element such that
%       the distance between a sample and any of its neighbors should be
%       less than thres.
%
%   [nbs, nbdists] = findnbs( ... );
%       also returns the distances from each sample to its neighbors.
%

%% parse and verify input arguments

assert(isnumeric(dists) && ndims(dists) == 2 && size(dists,1) == size(dists,2), ...
    'findnbs:invalidarg', ...
    'dists should be a square matrix.');

n = size(dists, 1);

assert(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1 && K <= n-1, ...
    'findnbs:invalidarg', ...
    'K should be an integer scalar with in [1, n-1]');
    
if nargin < 3
    thres = inf;
else
    assert(isnumeric(thres) && isscalar(thres) && isreal(thres) && thres > 0, ...
        'findobs:invalidarg', ...
        'thres should be a positive real value.');
end

out_dist = nargout >= 2;


%% main

% prevent from selecting self

dists(1:(n+1):n*n) = inf;

% sort and select

[sdists, si] = sort(dists, 1);
sdists = sdists(1:K, :);
si = si(1:K, :);

% output as cell array

nbs = cell(n, 1);

if ~out_dist
    for i = 1 : n
        csi = si(:, i);
        csd = sdists(:, i);
        
        nbs{i} = csi(csd < thres);
    end
else
    nbdists = cell(n, 1);
    
    for i = 1 : n
        csi = si(:, i);
        csd = sdists(:, i);
        
        nbs{i} = csi(csd < thres);
        nbdists{i} = csd(csd < thres);        
    end    
end
