function g = make_gr(dty, n, s, t, op)
% Makes a graph with given edges
%
%   g = make_gr('d', n, s, t);
%
%       makes a directed graph with given edges.
%
%       Here, n is the number of nodes, s is the source nodes, t is
%       the target nodes. s and t should be vectors of the same size.
%
%   g = make_gr('u', n, s, t);
%
%       makes an undirected graph with given edges.
%
%       For each pair of nodes, only one edge needs to be input.
%
%   g = make_gr('d', n, s, t, 'nbs');
%   g = make_gr('u', n, s, t, 'nbs');
%
%       makes a graph with neighborhood system established.
%
%   In the output, g is a struct with the following fields:
%
%   - tag:      a char string: 'gr'
%   - dty:      the direction of type: 'd' or 'u'
%   - n:        the number of nodes
%   - m:        the number of edges
%   - edges:    the list of edges (as columns): 2 x m' matrix 
%   - has_nbs:  whether neighborhood system is set up.
%
%   If has_nbs is true, it also has the following fields:
%   
%   - o_nbs:    the list of outgoing neighbors: 1 x m'
%   - o_eds:    the list of outgoing edges:     1 x m'
%   - o_degs:   the outgoing degrees:           1 x n
%   - o_os:     the offsets of sections:        1 x n
%
%   Here, m' = m when dty == 'd' and m' == 2*m when dty == 'u'
%

% Created by Dahua Lin, on Oct 28, 2011
%

%% verify input arguments

if ~(ischar(dty) && isscalar(dty) && (dty == 'd' || dty == 'u'))
    error('make_gr:invalidarg', 'dty should be either ''d'' or ''u''.');
end

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
    error('make_gr:invalidarg', 'n should be a non-negative integer number.');
end

if ~(isnumeric(s) && isvector(s) && isnumeric(t) && isvector(t) && ...
    isequal(size(s), size(t)))
    error('make_gr:invalidarg', ...
        's and t should be a numeric vectors of the same size.');
end

if size(s, 1) > 1; s = s.'; end  % turn to row vector
if size(t, 1) > 1; t = t.'; end

if nargin < 5
    use_nbs = false;
else
    if ~(ischar(op) && strcmpi(op, 'nbs'))
        error('make_gr:invalidarg', 'The 5th argument is invalid.');
    end
    use_nbs = true;
end

%% main

m = numel(s);

g.tag = 'gr';
g.dty = dty;
g.n = int32(n);
g.m = int32(m);

s = int32(s);
t = int32(t);

if dty == 'd'
    g.edges = [s; t];
else
    g.edges = [s t; t s];
end
g.has_nbs = false;

% setup neighborhood

if use_nbs
    g = gr_nbs(g);
end


