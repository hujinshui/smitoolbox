function [Gs, ws, Ws] = gr_sym(G, w, option)
% Make an undirected graph from a directed one via symmetrization
%
%   Gs = gr_sym(G, [], 'max');
%   Gs = gr_sym(G, [], 'min');
%
%       Make an undirected graph Gs from a directed graph G through
%       symmetrization. If the 2nd argument is 'max', then the edge
%       {i, j} will be preserved when either (i, j) or (j, i) is in G.
%       If the 2nd argument is 'min', then it will be preserved if
%       both (i, j) and (j, i) are in G.
%
%       G can be either a graph struct or a degree-bounded graph struct.
%
%   [Gs, ws] = gr_sym(G, w, 'avg');
%   [Gs, ws] = gr_sym(G, w, 'min');
%   [Gs, ws] = gr_sym(G, w, 'max');
%
%       Make an undirected weighted graph Gs from a directed weighted 
%       graph G through symmetrization, according to the 2nd argument.
%
%       The meaning of the 2nd argument
%       'avg':      Gs(i,j) = G(i,j) + G(j,i) / 2
%       'max':      Gs(i,j) = max(G(i,j), G(j,i))
%       'min':      Gs(i,j) = min(G(i,j), G(j,i))
%
%   [Gs, ws, Ws] = gr_sym( ... );
%
%       additionally returns the corresponding adjacency weight matrix Ws.
%
%   Note that the number of edges can change, if the some edge weights
%   change from zero to non-zero, or vice versa.
%
%   The output graph has no neighborhood system, and one can call gr_nbs
%   to build it.
%

% Created by Dahua Lin, on Nov 30, 2011
%

%% verify input arguments

if is_gr(G) && G.dty == 'd'
    if ~isempty(w)
        if ~(isnumeric(w) && isvector(w) && numel(w) == G.m)
            error('gr_sym:invalidarg', ...
                'w should be a numeric vector of length m.');
        end
    end
    
elseif is_gr_bnd(G)        
    if ~isempty(w)
        if ~(isnumeric(w) && isequal(size(w), [G.K, G.n]))
            error('gr_sym:invalidarg', ...
                'w should be a numeric matrix of size K x n');
        end
    end
    
else
    error('gr_sym:invalidarg', ...
        'G should be a directed graph struct or a degree-bounded graph struct.');
end


if ~(ischar(option))
    error('gr_sym:invalidarg', 'The 2nd argument is invalid.');
end

%% main

if isempty(w)    
    switch lower(option)
        case 'max'
            sym_fun = @or;
        case 'min'
            sym_fun = @and;
        otherwise
            error('gr_sym:invalidarg', 'The 2nd argument is invalid.');
    end
    
    W = gr_wmat(G);
    Ws = sym_fun(W, W.');
    [s, t] = find(Ws);
    i = find(s <= t);
    Gs = make_gr('u', G.n, s(i), t(i));
    ws = [];
    
else    
    switch lower(option)
        case 'max'
            sym_fun = @max;
        case 'min'
            sym_fun = @min;
        case 'avg'
            sym_fun = @(x,y) (x + y) / 2;
        otherwise
            error('gr_sym:invalidarg', 'The 2nd argument is invalid.');
    end

    W = gr_wmat(G, w);
    Ws = sym_fun(W, W.');
    [s, t, ws] = find(Ws);
    
    i = find(s <= t);
    s = s(i);
    t = t(i);
    ws = ws(i);
    Gs = make_gr('u', G.n, s, t);
    
end

