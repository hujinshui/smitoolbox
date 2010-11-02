function [s, t, w] = gridgraph1d(n, kernel, offsets)
% Construct locally connected graph based on 1D graph
%
%   [s,t,w] = gridgraph1d(n, kernel);
%       returns the edges of the graph among n nodes along a 1D grid.
%
%       kernel is a vector specifying the weights between pairs of
%       nodes with different distance. 
%
%       In the output, s, t, and w are respectively the sources, targets,
%       and weights of the undirected edges. Note that each edge appears
%       only once with s < t. Particularly, each vertex v is connected
%       to v-k and v+k with edge weight kernel(k).
%       
%   W = gridgraph1d(n, kernel, offsets)
%       returns the edges of the graph among n nodes along a 1D grid.
%
%       In this construction, we connect each vertex v to v+offsets(k)
%       and v-offsets(k) with edge weight kernel(k).
%
%       Note that gridgraph1d(n, kernel) is equivalent to
%       gridgraph1d(n, kernel, 1:h) where h is length(kernel).
%
%   Remarks
%   -------
%       With the resultant s, t, and w, one can make a graph struct
%       by invoking gr_edgelist or gr_adjlist.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 16, 2010
%       - Modified by Dahua Lin, on Sep 21, 2010
%       - Modified by Dahua Lin, on Nov 1, 2010
%

%% verify input arguments

if ~(isscalar(n) && isnumeric(n) && n == fix(n) && n >= 1)
    error('gridgraph1d:invalidarg', 'n should be a positive integer scalar.');
end

if ~(isfloat(kernel) && isreal(kernel) && isvector(kernel))
    error('gridgraph1d:invalidarg', 'kernel should be a real vector.');
end
h = numel(kernel);

if nargin < 3
    offsets = 1 : h;
else
    if ~(isnumeric(offsets) && isequal(size(offsets), size(kernel)) && ...
            all(offsets == fix(offsets)))
        error('gridgraph1d:invalidarg', ...
            'offsets should be an array of integers with same size of kernel.');
    end
end


%% main

kz = kernel == 0;
if any(kz)
    offsets = offsets(~kz);
    kernel = kernel(~kz);
    h = numel(kernel);
end

if h == 1
    s = 1 : n-offsets;
    t = 1+offsets : n;
    w = kernel(ones(1, n-offsets));
        
else
    ns = n * h - sum(offsets);
    s = zeros(1, ns);
    t = zeros(1, ns);
    w = zeros(1, ns);
    
    i = 0;
    for k = 1 : h
        o = offsets(k);
        
        si = i+1 : i+(n-o);        
        s(si) = 1 : n-o;
        t(si) = 1+o : n;
        w(si) = kernel(k);        
        i = i + (n-o);
    end
end


