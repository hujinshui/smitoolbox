function G = gr_adjlist(varargin)
% Construct or verify an adjacency list representation of a graph
%
%   In the representation used in smitoolbox, the adjacency list
%   representation is just a refined version of edgelist, such that
%   all edges are sorted according to its source vertex. And indexing
%   information is added to facilitate efficient retrieval of 
%   adjacent vertices for each vertex.
%
%   G = gr_adjlist(G);
%       if G is a valid adjacency list, it simply returns G, or if
%       G is an ordinary edge list, it upgrades it to an adjacency list.
%
%   G = gr_adjlist(A);
%       constructs an adjacency list from the adjacency matrix A.
%
%   G = gr_adjlist(n, [s, t]);
%   G = gr_adjlist(n, [s, t, w]);
%   G = gr_adjlist(n, s, t);
%   G = gr_adjlist(n, s, t, w);
%       constructs an adjacency list struct from given edges. 
%
%       In the input, n is the number of vertices, and s, t, w (optional)
%       are vectors of length m. For any i = 1, ..., m, s(i), t(i), and
%       w(i) are respectively the source vertex, target vertex, and weight
%       of the i-th object.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 9, 2010
%


%% main

if nargin == 1 && isstruct(varargin{1})
    
    G = varargin{1};
    
    if isfield(G, 'tag')
        
        if strcmp(G.tag, 'gr_edgelist')
            G = gr_adjlist_cimp(G);        
        elseif strcmp(G.tag, 'gr_adjlist') || strcmp(G.tag, 'gr_bdir_adjlist')
            return;
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
    
else
    G = gr_edgelist(varargin{:});    
    G = gr_adjlist_cimp(G);
    
end



