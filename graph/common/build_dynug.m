function [G, Gc] = build_dynug(nfrms, int_graphs, cross_edges)
% Builds a dynamic undirected graph across multiple frames
%
%   G = build_dynug(nfrms, int_graphs, cross_edges);
%
%       The function builds a graph over multiple frames. In this model,
%       each frame has its own internal graph, and there are edges
%       connecting the vertices in consecutive frames.
%
%       Input arguments:
%
%       - nfrms:        The number of frames
%
%       - int_graphs:   the internal graph of each frame. It can be 
%                       an undirected gr_edgelist object (it all frames
%                       have the graphical structure internally),
%                       or a cell array of such objects. Note that
%                       different frames need not have the same number
%                       of vertices.
%
%       - cross_edges:  the edges connecting between frames.
%                       If the cross-frame connection structure is the 
%                       same for every pair of consecutive frames,
%                       this can be a matrix with three columns. If not,
%                       it can be input as a cell array of such matrices.
%
%       The output is an undirected gr_adjlist object, which is over
%       all frames. The vertices are re-indexed in natural order.
%       The edges are re-indexed, such that the within-frame edges
%       come first, and the cross-frame edges come next.
%
%   {G, Gc] = build_dynug( ... );
%       
%       additionally returns Gc, which is a directed graph that contains 
%       only cross-frame edges.
%

%   Created by Dahua Lin, on Feb 16, 2011
%

%% verify input arguments

if ~(isnumeric(nfrms) && isscalar(nfrms) && nfrms == fix(nfrms) && nfrms >= 1)
    error('build_dynug:invalidarg', 'nfrms should be a positive integer scalar.');
end
nfrms = double(nfrms);

% int_graphs

if isa(int_graphs, 'gr_edgelist')
    
    if int_graphs.dtype ~= 'u'
        error('build_dynug:invalidarg', 'int_graphs should be undirected.');
    end
    shared_int = true;
    
    nvs = constmat(1, nfrms, int_graphs.nv);
    
elseif iscell(int_graphs)
    
    if numel(int_graphs) ~= nfrms
        error('build_dynug:invalidarg', ...
            'The number of cells in int_graphs is not equal nfrms.');
    end
    shared_int = false;
    
    nvs = zeros(1, nfrms);
    for k = 1 : nfrms
        g = int_graphs{k};
        
        if ~(isa(g, 'gr_edgelist') && g.dtype == 'u')
            error('build_dynug:invalidarg', ...
                'Each graph in int_graphs should be undirected.');
        end
        nvs(k) = g.nv;                
    end
    
else
    error('build_dynug:invalidarg', 'int_graphs is invalid.');
end

% cross_edges

if isnumeric(cross_edges)
    
    if ~(ndims(cross_edges) == 2 && size(cross_edges, 2) == 3)
        error('build_dynug:invalidarg', ...
            'The size of cross_edges matrix is invalid.');
    end
    shared_cro = true;
    
elseif iscell(cross_edges)
    
    if numel(cross_edges) ~= nfrms - 1
        error('build_dynug:invalidarg', ...
            'The number of cells in cross_edges is not equal (nfrms-1).');
    end
    
    for k = 1 : nfrms-1
        E = cross_edges{k};
        
        if ~(isnumeric(E) && ndims(E) == 2 && size(E,2) == 3)
            error('build_dynug:invalidarg', ...
                'Each cell in cross_edges should be an ne x 3 matrix.');
        end
    end
    shared_cro = false;
end

%% main

% re-index vertices

vbs = [0, cumsum(nvs(1:end-1))];

% collect internal edges

ess_int = cell(nfrms, 1);
ets_int = cell(nfrms, 1);
ews_int = cell(nfrms, 1);

for k = 1 : nfrms
    
    if shared_int
        g = int_graphs;
    else
        g = int_graphs{k};
    end
    
    ess_int{k} = g.es(1:g.ne) + (1 + vbs(k));
    ets_int{k} = g.et(1:g.ne) + (1 + vbs(k));
    ews_int{k} = g.ew(1:g.ne);
end

% collect cross-frame edges

ess_cro = cell(nfrms-1, 1);
ets_cro = cell(nfrms-1, 1);
ews_cro = cell(nfrms-1, 1);

for k = 1 : nfrms-1
    
    if shared_cro
        E = cross_edges;
    else
        E = cross_edges{k};
    end
    
    ess_cro{k} = E(:,1) + vbs(k);
    ets_cro{k} = E(:,2) + vbs(k+1);
    ews_cro{k} = E(:,3);
end
    
% make graph

tot_nv = sum(nvs);

s_a = vertcat(ess_int{:}, ess_cro{:});
t_a = vertcat(ets_int{:}, ets_cro{:});
w_a = vertcat(ews_int{:}, ews_cro{:});

G = gr_adjlist.from_edges('u', tot_nv, s_a, t_a, w_a);

if nargout >= 2
    
    s_c = vertcat(ess_cro{:});
    t_c = vertcat(ets_cro{:});
    w_c = vertcat(ews_cro{:});
    
    Gc = gr_adjlist.from_edges('d', tot_nv, s_c, t_c, w_c);
end


