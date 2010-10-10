function A = edgelist_to_adjmat(G, op)
% Converts an edge list to an adjacency matrix
%
%   A = edgelist_to_adjmat(G);
%   A = edgelist_to_adjmat(G, 'full');
%
%       G is a struct to represent an edge list or an adjacency list.
%       
%       In the output, A is an n x n matrix, such that only the elements
%       at A(G.s(i), G.t(i)) are non-zero. Moreover, if G.w is non-empty,
%       then A is a numeric matrix with the same type as G.w, and it 
%       has A(G.s(i), G.t(i)) = G.w(i) for each i = 1, ..., m.
%
%       By default, a sparse matrix is created. If ones want to create
%       a full matrix, then set the 2nd argument to 'full'.
%

% Created by Dahua Lin, on Oct 9, 2010


%% verify

if ~(numel(G) == 1 && isstruct(G))
    error('edgelist_to_adjmat:invalidarg', ...
        'G should be a edge list struct that represents an edge list.');
end

if nargin >= 2
    if ~strcmp(op, 'full')
        error('edgelist_to_adjmat:invalidarg', ...
            'The 2nd argument can only be ''full''.');
    end
    use_full = 1;
else
    use_full = 0;
end
    
    
%% main

n = double(G.n);

if use_full

    s = G.s;
    t = G.t;
    
    if isempty(G.w)
        A = false(n, n);
        A(s + t * n + 1) = 1;
    else
        A = zeros(n, n, class(G.w));
        A(s + t * n + 1) = G.w;
    end
        
else
    
    s = double(G.s+1);
    t = double(G.t+1);
    
    if isempty(G.w)
        A = sparse(s, t, true, n, n);
    else
        A = sparse(s, t, double(G.w), n, n);
    end    
    
end



