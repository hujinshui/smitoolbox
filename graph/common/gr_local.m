function G = gr_local(siz, rgn)
% Create a degree-bounded graph with links between local neighbors
%
%   G = gr_local(n, r);
%
%       Creates a locally connected graph over one-dimensional
%       grid. Here, n is the number of nodes, and r is the local range.
%       In this graph, each node is connected to a neighboring node,
%       within distance r.
%   
%   G = gr_local([m, n], r);
%   G = gr_local([m, n], [yr, xr]);
%
%       Creates a locally connected graph over two-dimensional
%       grid. Here, the grid size is m rows and n columns, and yr and xr
%       are respectively the local range along y and x direction.
%
%       If xr == yr, one can simply input a single value r that equals
%       xr and yr.
%
%   G = gr_local([m, n], ker);
%
%       Creates a locally connected graph over two-dimensional
%       grid. Here, ker is a matrix specifying the connection pattern, 
%       the node (i, j) is connected to (i+di, j+dj) if ker(di', dj') 
%       is non-zero. Here, ker is a matrix of size [2*yr+1, 2*xr+1], and
%       di' = di+(yr+1) and dj' = dj+(dx+1).
%

% Created by Dahua Lin, on Nov 30, 2011
%

%% verify input arguments

if ~(isnumeric(siz) && ndims(siz) == 2 && all(siz == fix(siz)) )
    error('gr_local:invalidarg', 'The 1st argument is invalid.');
end

d = numel(siz);
if ~(d == 1 || d == 2)
    error('gr_local:invalidarg', 'Only 1D or 2D grid is allowed.');
end

if d == 1
    n = siz;
    r = rgn;
    if ~(isnumeric(r) && isscalar(r) && r == fix(r) && r >= 0)
        error('gr_local:invalidarg', 'r should be a positive integer.');
    end
    
else % d == 2
    m = siz(1);
    n = siz(2);
    
    if ndims(rgn) ~= 2
        error('gr_local:invalidarg', 'The 2nd argument is invalid.');
    end
    
    if isnumeric(rgn) && numel(rgn) <= 2
        if isscalar(rgn)
            yr = rgn;
            xr = rgn;
        else
            yr = rgn(1);
            xr = rgn(2);
        end
        ker = true(2*yr+1, 2*xr+1);
        ker(yr+1, xr+1) = 0;
        
    elseif islogical(rgn)        
        ker = rgn;
        [kh, kw] = size(ker);
        if ~(rem(kh, 2) == 1 && rem(kw, 2) == 1)
            error('gr_local:invalidarg', ...
                'The dimensions of the kernel should be odd numbers.');
        end
        yr = (kh - 1) / 2;
        xr = (kw - 1) / 2;
        
    else
        error('gr_local:invalidarg', 'The 2nd argument is invalid.');
    end    
end

%% main

if d == 1    
    delta = [-r:-1, 1:r].'; 
    nbs = bsxfun(@plus, delta, 1:n);
    nbs(nbs < 0 | nbs > n) = 0;
else
    [dx, dy] = meshgrid(-xr:xr, -yr:yr);
    [x, y] = meshgrid(1:n, 1:m);
    
    dx = dx(ker);
    dy = dy(ker);
    x = reshape(x, 1, m * n);
    y = reshape(y, 1, m * n);
    
    nx = bsxfun(@plus, dx, x);
    ny = bsxfun(@plus, dy, y);
    
    is_out = nx < 1 | nx > n | ny < 1 | ny > m;
    nbs = ny + (nx - 1) * m;
    nbs(is_out) = 0;
end
    
G = make_gr_bnd(nbs);
    
