function W = mrf2d(siz, kernel)
% Construct a second-order MRF between nodes with 2D index set
%
%   W = mrf2d([m, n], kernel);
%       returns the affinity matrix of the constructed MRF among
%       m x n nodes with 2D index set, here m is the number of rows
%       and n is the number of columns.
%
%       kernel should be a (h+1) x (w+1) matrix, where h < m and w < n.
%       Then the weights between the node at (i, j) and the node at
%       (i+k, j+l) is kernel(k+1, l+1), when k <= h and l <= w.
%       kernel(1, 1) should be 0, such that no link is added between
%       a node and itself.
%
%       The weights between other pairs of nodes are set to zeros.       
%
%       In the output, W is a sparse matrix of size N x N, where 
%       N = m x n is the number of nodes.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 16, 2010
%       - Modified by Dahua Lin, on Sep 21, 2010
%           - use new implementation
%

%% verify input arguments

if ~(isnumeric(siz) && numel(siz) == 2 && all(siz == fix(siz) & siz >= 1))
    error('mrf2d:invalidarg', ...
        'The first argument [m, n] should be a pair of positive integer scalars.');
end

h = siz(1);
w = siz(2);
N = h * w;

if ~(isfloat(kernel) && isreal(kernel) && ndims(kernel) == 2)
    error('mrf2d:invalidarg', 'kernel should be a real matrix.');
end

if kernel(1) ~= 0
    error('mrf2d:invalidarg', 'kernel(1, 1) must be zero.');
end

%% main

[kh, kw] = size(kernel);

dX = repmat(0:kw-1, kh, 1);
dY = repmat((0:kh-1)', 1, kw);

nk = numel(kernel);

I = cell(1, nk);
J = cell(1, nk);
V = cell(1, nk);

for k = 1 : nk
    
    kv = kernel(k);
    if kv == 0
        continue;
    end
  
    dx = dX(k);
    dy = dY(k);
    
    x0 = repmat(1:w-dx, h-dy, 1);
    if dx ~= 0 
        x1 = repmat(1+dx:w, h-dy, 1);
    end
        
    y0 = repmat((1:h-dy).', 1, w-dx);
    if dy ~= 0
        y1 = repmat((1+dy:h).', 1, w-dx);
    end
    
    if dx == 0 
        [i0, j0] = make_links(w, h, x0, x0, y0, y1);
        [i1, j1] = make_links(w, h, x0, x0, y1, y0);
        
        I{k} = [i0; i1];
        J{k} = [j0; j1];
        
    elseif dy == 0
        [i0, j0] = make_links(w, h, x0, x1, y0, y0);
        [i1, j1] = make_links(w, h, x1, x0, y0, y0);
        
        I{k} = [i0; i1];
        J{k} = [j0; j1];
        
    else
        [i0, j0] = make_links(w, h, x0, x1, y0, y1);
        [i1, j1] = make_links(w, h, x0, x1, y1, y0);
        [i2, j2] = make_links(w, h, x1, x0, y0, y1);
        [i3, j3] = make_links(w, h, x1, x0, y1, y0);
        
        I{k} = [i0; i1; i2; i3];
        J{k} = [j0; j1; j2; j3];
    end

    V{k} = constmat(numel(I{k}), 1, kv);
end
    
I = vertcat(I{:});
J = vertcat(J{:});
V = vertcat(V{:});

W = sparse(I, J, V, N, N);  


%% subfunctions

function [i, j] = make_links(w, h, xi, xj, yi, yj) %#ok<INUSL>

i = yi(:) + h * (xi(:) - 1);
j = yj(:) + h * (xj(:) - 1);


