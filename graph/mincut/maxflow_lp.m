function [fg, fs, ft] = maxflow_lp(g, sc, tc, solver)
% Solve the maximum flow problem using LP solver
%
%   [fg, fs, ft] = maxflow_lp(g, sc, tc);
%   [fg, fs, ft] = maxflow_lp(g, sc, tc, solver);
%
%       solves a maximum flow problem using LP solver.
%
%       Input arguments:
%       - g:    the capacity graph [gr_edgelist or gr_adjlist with n nodes]
%       - sc:   the capacity from source to nodes [a vector of length n]
%       - tc:   the capacity from nodes to target [a vector of length n]
%       - solver: the solver function handle to solve the LP problem.
%                 If omitted, this function uses @mstd_solve by default.
%
%       Output arguments:
%       - fg:   the flow along edges in g [m x 1 vector]
%       - fs:   the flow from source to nodes [n x 1 vector]
%       - ft:   the flow from nodes to target [n x 1 vector]
%

% Created by Dahua Lin, on Apr 17, 2011
%

%% verify input arguments

% The maxflow_lp_prob invoked below will check g, sc, and tc.

if nargin < 4
    solver = @mstd_solve;
else
    if ~isa(solver, 'function_handle')
        error('maxflow_lp:invalidarg', 'solver must be a function handle.');
    end
end

%% main

% construct problem

P = maxflow_lp_prob(g, sc, tc);

n = g.nv;
m = numel(g.es);

si = find(sc > 0);
ti = find(tc > 0);

m0 = numel(si);
m1 = numel(ti);

if P.d ~= m + m0 + m1
    error('maxflow_lp:rterror', 'Inconsistent solution dimension.');
end

% solve

x = solver(P);

% extract results

fg = x(1:m);
fs = zeros(n, 1);
fs(si) = x(m+1 : m+m0);
ft = zeros(n, 1);
ft(ti) = x(m+m0+1 : m+m0+m1);


% post-process the results
% to ensure that for each pair of edges, at least one has zero flow

% find out the flow of reverse edge

if ~g.is_directed % undirected    
    ne = g.ne;
    fr = [fg(ne+1 : 2*ne); fg(1 : ne)]; 
    
else  % directed
    sv = double(g.source_vs);
    tv = double(g.target_vs);
    
    F = sparse(sv, tv, fg, n, n);
    fr = full(F(tv + n * (sv-1)));
end
    
fg = fg - min(fg, fr);





