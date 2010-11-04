function A = edgelist_to_adjmat(G, varargin)
% Converts an edge list to an adjacency matrix
%
%   A = edgelist_to_adjmat(G);
%   A = edgelist_to_adjmat(G, 'full');
%   A = edgelist_to_adjmat(G, dty);
%   A = edgelist_to_adjmat(G, dty, 'full');
%
%       G is a struct to represent an edge list or an adjacency list.
%       
%       if dty is not specified, by default it uses G.dty if it is
%       available. If G.dty is not available, then it uses 'd' by default.
%
%       The function creates an n x n matrix A, with
%       A(G.s(i), G.t(i)) = G.w(i) (for weighted graph) or true
%       (for unweighed graph). If dty is 'u', then 
%       A(G.t(i), G.s(i)) are also set to the same value.
%       All other entries in A are set to zero.       
%
%       By default, a sparse matrix is created. If ones want to create
%       a full matrix, then set the 2nd argument to 'full'.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 9, 2010
%       - Modified by Dahua Lin, on Oct 31, 2010
%           - support undirected graph


%% verify

if ~(numel(G) == 1 && isstruct(G))
    error('edgelist_to_adjmat:invalidarg', ...
        'G should be a edge list struct that represents an edge list.');
end


use_full = 0;
dty = [];

if ~isempty(varargin)
    nv = numel(varargin);
    for i = 1 : nv
        switch varargin{i}
            case 'full'
                use_full = 1;
            case 'd'
                dty = 'd';
            case 'u'
                dty = 'u';
        end
    end
end
                    
add_re = ~(isfield(G, 'dty') && isequal(G.dty, 'u')) && isequal(dty, 'u');
    
%% main

n = double(G.n);

if use_full

    s = G.s;
    t = G.t;
    
    if isempty(G.w)
        A = false(n, n);
        A(s + t * n + 1) = 1;
        if add_re
            A(t + s * n + 1) = 1;
        end        
    else
        A = zeros(n, n, class(G.w));
        A(s + t * n + 1) = G.w;
        if add_re
            A(t + s * n + 1) = G.w;
        end
    end
        
else
    
    s = double(G.s+1);
    t = double(G.t+1);
    
    if ~add_re
        if isempty(G.w)
            A = sparse(s, t, true, n, n);
        else
            A = sparse(s, t, double(G.w), n, n);
        end    
    else
        if isempty(G.w)
            A = sparse([s; t], [t; s], true, n, n);
        else
            A = sparse([s; t], [t; s], double([G.w; G.w]), n, n);
        end 
    end
end



