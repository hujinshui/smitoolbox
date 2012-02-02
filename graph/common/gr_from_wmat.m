function [g, w] = gr_from_wmat(W, dt, op)
% Construct a graph struct from an edge-weight matrix
%
%   [g, w] = gr_from_wmat(W, 'u');
%   [g, w] = gr_from_wmat(W, 'd');
%   [g, w] = gr_from_wmat(W, 'u', 'nbs');
%   [g, w] = gr_from_wmat(W, 'd', 'nbs');
%
%       Constructs a graph g from W, the edge weight matrix. 
%       
%       The second argument specifies whether to construct an directed
%       or undirected graph, depending on its value ('d' or 'u').
%       For undirected graph, W needs to be a symmetric matrix.
%
%       If the third argument is given as 'nbs', the neighborhood struct
%       will also be established.
%
%       In the output, g is the graph struct, and w is the vector of 
%       edge weights.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 28, 2011
%

%% verify input arguments

n = size(W, 1);
if ~((isnumeric(W) || islogical(W)) && ndims(W) == 2 && size(W,2) == n)
    error('gr_from_wmat:invalidarg', ...
        'W should be a square matrix (numeric or logical).');
end

if ~(ischar(dt) && (dt == 'u' || dt == 'd'))
    error('gr_from_wmat:invalidarg', ...
        'dt should be a char which equals ''u'' or ''d''.');
end

if nargin >= 3
    if ~(ischar(op) && strcmpi(op, 'nbs'))
        error('gr_from_wmat:invalidarg', ...
            'The 3rd argument should be ''nbs'' if given.');
    end
    use_nbs = true;
else
    use_nbs = false;
end


%% main

[s, t, w] = find(W);

if dt == 'u'
    r = find(s <= t);
    s = s(r);
    t = t(r);
    w = w(r);
end

if use_nbs
    g = make_gr(dt, n, s, t, 'nbs');
else
    g = make_gr(dt, n, s, t);
end





