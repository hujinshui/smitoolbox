function [labels, maxflow] = graphcut(n, edges, tweights, nweights)
%GRAPHCUT Performs Boycov's graphcut on a given weighted graph
%
%   [labels, maxflow] = graphcut(n, edges, tweights, nweights)
%
%   Input:
%       - n:        the number of nodes
%       - edges:    the array of edges, 
%                   It is a 2 x m array. m is the number of edges.
%                   Each column of edges is [isource; isink], which are indices
%                   of end nodes of an edge.
%       - tweights: the terminal weights, which equals to the difference between
%                   the weight of being tied to source and that of being tied
%                   to sink.
%       - nweights: the neighboring weights. 
%                   
%   Output: 
%       - labels:   the labels of each node given by the solution.
%                   It is a logical array of size 1 x n. 
%                   0: assigned to source part, 1: assigned to sink part.
%       - maxflow:  the value of maximum flow (minimum cut).
%
%   Remarks:
%       1. All edges are considered as symmetric. If you specify an edge (i, j)
%          with neighbor weight w, then essentially two edges will be added to
%          the graph: (i, j) and (j, i) of weight w.
%

%   Created by Dahua Lin, on Oct 11, 2009
%   Modified by Dahua Lin, on Nov 6, 2009 
%       - change weight type from int32 to double
%

%% verify input arguments

error(nargchk(4, 4, nargin));

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n > 0)
    error('graphcut:invalidarg', 'n should be a positive integer scalar.');
end

n = double(n);

if ~(isnumeric(edges) && ndims(edges) == 2 && size(edges, 1) == 2)
    error('graphcut:invalidarg', 'edges should be a 2 x m matrix.');
end
m = size(edges, 2);

if ~isa(edges, 'int32')
    edges = int32(edges);
end

if ~(isnumeric(tweights) && isequal(size(tweights), [1, n]))
    error('graphcut:invalidarg', 'tweights should be a numeric vector of size 1 x n.');
end

if ~isa(tweights, 'double')
    tweights = double(tweights);
end

if ~(isnumeric(nweights) && isequal(size(nweights), [1, m]))
    error('graphcut:invalidarg', 'nweights should be a numeric vector of size 1 x m.');
end

if ~isa(nweights, 'double')
    nweights = double(nweights);
end


%% main

[labels, maxflow] = graphcut_cimp(n, edges, tweights, nweights);



