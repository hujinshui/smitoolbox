function G = gr_adjlist(varargin)
% Construct or verify an adjacency list representation of a graph
%
%   In the representation used in smitoolbox, the adjacency list
%   representation is an augmented version of edgelist, with additional
%   data structure to represent neighborhood of each vertex.
%
%   G = gr_adjlist(G, dty);
%       verifies the validity of an adjacency list struct of specified
%       type, where dty can be 'd' or 'u', respectively indicating
%       directed or undirected graph.
%
%   G = gr_adjlist(edlist, dty);
%       constructs a adjacency list from an input edge list.
%       Here, dty can be 'd' or 'u', indicating whether the graph is
%       a directed or undirected graph. If dty is omitted, it uses 'd'.
%
%   G = gr_adjlist(A, dty);
%       constructs an adjacency list from the adjacency matrix A.
%       Here, dty can be 'd' or 'u', indicating whether the graph is
%       a directed or undirected graph. If dty is omitted, it uses 'd'.
%
%   G = gr_adjlist(dty, n, [s, t]);
%   G = gr_adjlist(dty, n, [s, t, w]);
%   G = gr_adjlist(dty, n, s, t);
%   G = gr_adjlist(dty, n, s, t, w);
%       constructs an adjacency list struct from given edges. 
%       Here, dty can be 'd' or 'u', indicating whether the graph is
%       a directed or undirected graph.
%
%       In the input, n is the number of vertices, and s, t, w (optional)
%       are vectors of length m. For any i = 1, ..., m, s(i), t(i), and
%       w(i) are respectively the source vertex, target vertex, and weight
%       of the i-th object.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 9, 2010
%       - Modified by Dahua Lin, on Oct 31, 2010
%           - the new struct preserves the input order of edges
%


%% main

if nargin <= 2 
    
    if isstruct(varargin{1})
    
        G = varargin{1};  
        if nargin < 2
            dty = [];
        end
    
        if isfield(G, 'tag')
            if strcmp(G.tag, 'gr_adjlist')
                if isempty(dty)
                    return;
                else
                    if ~isequal(G.dty, dty)
                        error('gr_adjlist:invalidarg', ...
                            'The direction type is not as expected.');
                    end
                end
                
                
            elseif strcmp(G.tag, 'gr_edgelist')
                if isempty(dty)
                    dty = 'd';
                end
                G = make_adjlist(G, dty);
                
            else
                G = [];
            end
        else
            G = [];
        end
        
        if isempty(G)
            error('gr_adjlist:invalidarg', ...
                'G is not a valid struct for graph representation.');
        end
    
    elseif isnumeric(varargin{1})
        
        A = varargin{1};
        if nargin < 2
            dty = 'd';
        else
            dty = varargin{2};
        end
                    
        G = gr_edgelist(A, dty);
        G = make_adjlist(G, dty);
    end                
    
else
        
    dty = varargin{1};
    
    G = gr_edgelist(varargin{2:end});    
    G = make_adjlist(G, dty);
    
end


%% The core function for creating adjlist struct

function G = make_adjlist(G, dty)

if isequal(dty, 'd')    
    [G.o_ds, G.o_os, G.o_es, G.o_ns] = build_nbs(G.n, G.s, G.t);
    
elseif isequal(dty, 'u')
    [G.o_ds, G.o_os, G.o_es, G.o_ns] = build_nbs(G.n, [G.s; G.t], [G.t; G.s]);
    G.o_es = rem(G.o_es, int32(G.m));
    
else
    error('gr_adjlist:invalidarg', ...
        'The direction type is invalid.');
end

G.dty = dty;
G.tag = 'gr_adjlist';


function [ds, os, es, ns] = build_nbs(n, s, t)

[s, es] = sort(s);
ns = t(es);

p = valueseg(s);
pds = [diff(p); numel(s)-p(end)+1];

ds = zeros(n, 1);
ds(s(p)+1) = pds;

os = [0; cumsum(ds(1:end-1))];

ds = int32(ds);
os = int32(os);
es = int32(es) - 1;
ns = int32(ns);

