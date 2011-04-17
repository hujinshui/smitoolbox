function P = maxflow_lp_prob(g, sc, tc)
% Constructs the LP maxflow problem
%
%   The maximum flow problem is formulated as
%
%       maximize sum_v c(s, v) * f(s, v);
%
%           s.t. f(s, v) + sum_u f(u, v) = f(v, t) +  sum_w f(v, w)
%
%                0 <= f(s, v) <= c(s, v)  for each v with c(s, v) > 0
%                0 <= f(v, t) <= c(v, t)  for each v with c(v, t) > 0
%                0 <= f(u, v) <= c(u, v)  for each (u, v) 
%       
%
%   P = maxflow_lp_prob(g, sc, tc);
%       
%       Constructs the primal LP problem for a maximum flow problem, 
%       as formulated above.
%
%       Input arguments:
%       - g:    the capacity graph [gr_edgelist or gr_adjlist with n nodes]
%       - sc:   the capacity from source to nodes [a vector of length n]
%       - tc:   the capacity from nodes to target [a vector of length n]
%
%       This functions returns a lp_problem struct. Suppose there are
%       m edges in g, and m0 source-node edges (non-zeros in sc), and
%       m1 node-target edges (non-zeros in tc).
%
%       Note that the constructed problem is a minimization problem,
%       whose objective coefficients have negative signs.
%
%       Let x be the solution, then its length is m+m0+m1, and
%       x(1:m) are the flows between nodes in g, 
%       x(m+1:m+m0) are the flows from source to nodes, and
%       x(m+m0+1:m+m0+m1) are the flows from nodes to target.
%

% Created by Dahua Lin, on Apr 17, 2011
%

%% verify input arguments

if ~( isa(g, 'gr_adjlist') && ~isempty(g.ew) )
    error('maxflow_lp_prob:invalidarg', ...
        'g should be an object of class gr_adjlist with weighted edges.');
end

n = g.nv;

if ~(isnumeric(sc) && isreal(sc) && isvector(sc) && length(sc) == n)
    error('maxflow_lp_prob:invalidarg', ...
        'sc should be a real vector of length n.');
end
if ~(isnumeric(tc) && isreal(tc) && isvector(tc) && length(tc) == n)
    error('maxflow_lp_prob:invalidarg', ...
        'tc should be a real vector of length n.');
end

if ~isa(sc, 'double'); sc = double(sc); end
if ~isa(tc, 'double'); tc = double(tc); end

if size(sc, 2) > 1; sc = sc.'; end
if size(tc, 2) > 1; tc = tc.'; end


%% main

% prepare

u = double(g.source_vs);
v = double(g.target_vs);
w = g.weights;
if ~isa(w, 'double'); w = double(w); end

si = find(sc > 0);
ti = find(tc > 0);

m = numel(u);
m0 = numel(si);
m1 = numel(ti);

% objective

f = - [zeros(m, 1); ones(m0, 1); zeros(m1, 1)];

% network constraints

M = sparse([u; v], [(1:m)'; (1:m)'], [-ones(m, 1); ones(m, 1)], n, m);
M0 = sparse(si, 1:m0, 1, n, m0);
M1 = sparse(ti, 1:m1, -1, n, m1);

Ae =[M, M0, M1];
be = zeros(n, 1);

% capacity constraints (upper bounds)

ub = [w; sc(si); tc(ti)];

% make problem struct

P = lp_problem(f, [], [], Ae, be, 0, ub);


