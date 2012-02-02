function [x, energy] = minenergy_graphcut(g, w, cp, cn)
%MINENERGY_GRAPHCUT Energy minimization via Graph-cut
%
%   [x, energy] = MINENERGY_GRAPHCUT(g, w, cp, cn);
%
%       Solves the following energy minimization problem using minimum
%       graph cut algorithm.
%
%           minimize E(L) = sum_{i in V} cp(i) * I(x_i == 1)
%                         + sum_{i in V} cn(i) * I(x_i == -1)
%                         + sum_{(i,j) in E} w_ij I(x_i <> x_j)
%
%           x_i can take value of either 1 or -1.
%           
%       Here, all coefficients cp(i), cn(i), and w_ij need to be
%       non-negative. 
%
%
%       Input arguments:
%       - g:        The underlying undirected graph
%       
%       - w:        The edge weights (for encouraging smoothness), which
%                   can be either a scalar or a vector of length g.m.
%
%       - cp:       The costs of assigning 1 to each vertex.
%                   (either a scalar or a vector of length g.n)
%
%       - cn:       The costs of assigning -1 to each vertex.
%                   (either a scalar or a vector of length g.n)
%
%       Output arguments:
%       - x:        The solution vector [1 x n int32].
%                   x(i) can be either 1, -1 or 0.
%                   x(i) == 0 indicates that one can assign either 1 or
%                   -1 to x(i).
%                   
%       - energy:   The minimum energy value.
%


% Created by Dahua Lin, on Jan 29, 2011
%

%% verify input arguments

if ~(is_gr(g) && isequal(g.dty, 'u'))
    error('minenergy_graphcut:invalidarg', ...
        'g should be a graph struct (for undirected graph).');
end

w = check_cvec('w', w, g.m);
cp = check_cvec('cp', cp, g.n);
cn = check_cvec('cn', cn, g.n);

%% main

[x, energy] = kolmogorov_mincut(g, w, cn, cp);


%% auxiliary functions

function v = check_cvec(name, v, n)

if ~(isfloat(v) && (isscalar(v) || (isvector(v) && isreal(v) && numel(v) == n)))
    error('minenergy_graphcut:invalidarg', ...
        '%s should be either a scalar or a vector of proper length.', name);
end

if ~all(v >= 0)
    error('minenergy_graphcut:invalidarg', ...
        'All values in %s need to be non-negative.', name);
end

if isscalar(v)
    v = v(ones(1, n));
end

