function G = gr_edgelist(varargin)
% Construct or verify an edge list representation of a graph
%
%   G = gr_edgelist(G);
%       verifies whether G is a valid edge list array. If not, it raises
%       an error.
%
%       A gr_edgelist struct contains the following fields:
%       - tag:  a string indicating the type of the graph
%       - n:    the number of vertices
%       - m:    the number of edges
%       - s:    the source vertices of all edges [m x 1 int32 zero-based]
%       - t:    the target vertices of all edges [m x 1 int32 zero-based]
%       - w:    the edge weights [empty or m x 1 numeric]
%
%   G = gr_edgelist(A);
%   G = gr_edgelist(A, 'u');
%       constructs an edge list representation from an adjacency matrix.
%       
%       A here can be either a logical or numeric array. If W is logical,
%       then it creates an edge list with unweighted edges, otherwise
%       it creates an edge list with weighted edges.
%
%       When the second argument is 'u', only those pairs (s, t) with
%       s < t are added.
%
%   G = gr_edgelist(n, [s, t]);
%   G = gr_edgelist(n, [s, t, w]);
%   G = gr_edgelist(n, s, t);
%   G = gr_edgelist(n, s, t, w);
%       constructs an edge list representation from explicitly given
%       edges. 
%
%       Here, s, t are source and target node indices, and w are
%       corresponding edge weights.
%

%% main

% determine input types

if nargin == 1
    if isstruct(varargin{1})
        intype = 'g';
    elseif isnumeric(varargin{1})
        intype = 'a';
        dty = 'd';
    else
        error('gr_edgelist:invalidarg', 'The input arguments are invalid.');
    end
    
elseif nargin == 2
    if ~isnumeric(varargin{1})
        error('gr_edgelist:invalidarg', 'The input arguments are invalid.');
    end
    
    if ischar(varargin{2})
        intype = 'a';
        dty = varargin{2};
    else
        intype = 'e';
    end
    
elseif nargin == 3 || nargin == 4
    intype = 'l';
else
    error('gr_edgelist:invalidarg', ...
        'The number of input arguments are invalid.');
end
    

% do construction

if intype == 'g'
    
    G = varargin{1};
    if ~( isfield(G, 'tag') && ( ...
            strcmp(G.tag, 'gr_edgelist') || ...
            strcmp(G.tag, 'gr_adjlist')) )
        
        error('gr_edgelist:invalidarg', ...
            'G is not a valid graph edgelist.');
    end
    
elseif intype == 'a'
    
    A = varargin{1};
    G = from_amat(A, dty);
    
elseif intype == 'e'        
    
    n = varargin{1};
    E = varargin{2};
    
    if ~(ndims(E) == 2 && isnumeric(E))
        error('gr_edgelist:invalidarg', ...
            'The 2nd argument should be a numeric matrix.');
    end
    
    m = size(E, 1);
    
    if size(E,2) == 2
        s = E(:,1);
        t = E(:,2);
        w = [];
        
    elseif size(E,2) == 3
        s = E(:,1);
        t = E(:,2);
        w = E(:,3);
        
    else
        error('gr_edgelist:invalidarg', ...
            'The 2nd argument should have 2 or 3 columns.');
    end
    
    G = from_stw(n, m, s, t, w);
    
else
    
    n = varargin{1};
    s = varargin{2};
    t = varargin{3};        
    
    if nargin < 4
        w = [];
    else
        w = varargin{4};
    end
    
    [m, s, t, w] = verify_stw(s, t, w);
    
    G = from_stw(n, m, s, t, w);    
    
end
    
    
    
%% construction functions

function G = from_amat(A, dty)
        
n = size(A, 1);
if ~(ndims(A) == 2 && n == size(A,2))
    error('gr_edgelist:invalidarg', 'A should be a square matrix.');
end

if isnumeric(A)
    [s, t, w] = find(A);
else
    [s, t] = find(A);
    w = [];
end

if dty == 'u'
    se = find(s < t);
    s = s(se);
    t = t(se);
    
    if ~isempty(w)
        w = w(se);
    end
end

m = numel(s);

G = from_stw(n, m, s, t, w);


function [m, s, t, w] = verify_stw(s, t, w)

if ~(isnumeric(s) && isreal(s) && isvector(s))
    error('gr_edgelist:invalidarg', 's should be a real vector.');
end

if ~(isnumeric(t) && isreal(t) && isvector(t))
    error('gr_edgelist:invalidarg', 't should be a real vector.');
end

m = numel(s);

if issparse(s); s = full(s); end
if issparse(t); t = full(t); end

if size(s,2) > 1; s = s.'; end
if size(t,2) > 1; t = t.'; end

if ~isempty(w)
    if ~(isnumeric(w) && isreal(w) && isvector(w) && numel(w) == m)
        error('gr_edgelist:invalidarg', 'w should be a real vector of length m.');
    end
    
    if issparse(w); w = full(w); end
    if size(w,2) > 1; w = w.'; end
end


function G = from_stw(n, m, s, t, w)

G = struct('tag', 'gr_edgelist', ...
    'n', n, 'm', m, 's', int32(s)-1, 't', int32(t)-1, 'w', w);




