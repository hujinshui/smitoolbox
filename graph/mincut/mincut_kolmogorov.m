function [L, cutv] = mincut_kolmogorov(G, src_ws, sink_ws)
% Solve minimum cut of undirected graph using Kolmogorov's implementation
%
%   [L, cutv] = mincut_kolmogorov(G, src_ws, sink_ws);
%       solves the minimum cut of directed graph using Kolmogorov's 
%       implementation.
%
%       Input:
%       - G:        the undirected graph (in a form convertible to gr_edgelist)
%       - src_ws:   the edge weights of each vertex to source
%       - sink_ws:  the edge weights of each vertex to sink
%
%       Output:
%       - L:        the cut result. L is a logical vector, L(i) == 0
%                   means the i-th vertex is assigned to the source 
%                   part, and L(i) == 1 means that i-th vertex is 
%                   assigned to the sink part.
%       - cutv:     the cut value of the solution.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 15, 2010
%

%% verify input arguments

G = gr_edgelist(G);

if isempty(G.w);
    error('mincut_kolmogorov:invalidarg', ...
        'G should be a weighted graph.');
end

n = G.n;

if ~(isvector(src_ws) && numel(src_ws) == n && isnumeric(src_ws))
    error('mincut_kolmogorov:invalidarg', ...
        'src_ws should be an 1 x n numeric vector.');
end

if ~(isvector(sink_ws) && numel(sink_ws) == n && isnumeric(sink_ws))
    error('mincut_kolmogorov:invalidarg', ...
        'sink_ws should be an 1 x n numeric vector.');
end

if issparse(src_ws);  src_ws = full(src_ws); end
if issparse(sink_ws); sink_ws = full(sink_ws); end

wc = class(G.w);
if ~isa(src_ws, wc); src_ws = cast(src_ws, wc); end
if ~isa(sink_ws, wc); sink_ws = cast(sink_ws, wc); end


%% main

[L, cutv] = mincut_kolmogorov_cimp(G, src_ws, sink_ws);

