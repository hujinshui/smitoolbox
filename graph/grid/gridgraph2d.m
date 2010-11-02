function [s,t,w] = gridgraph2d(siz, kernel)
% Constructs a locally connected graph based on a 2D grid
%
%   [s,t,w] = gridgraph2d([m, n], kernel);
%       constructs a graph based on a 2D grid with m rows and n columns.
%
%       In input, kernel can be a matrix of size kh x kw. Then in the
%       constructed graph, the vertex at (i, j) are connected to 
%       (i+di, j+dj) with |di| < kh and |dj| < kw. Particularly,
%       The edge connecting (i,j) to (i+di, j+dj) has a weight given
%       by kernel(di+1, dj+1).
%
%       In output, s, t, and w are respectively the sources, targets,
%       and weights of the undirected edges.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 16, 2010
%       - Modified by Dahua Lin, on Sep 21, 2010
%       - Modified by Dahua Lin, on Nov 1, 2010
%

%% verify input arguments

if ~(isnumeric(siz) && numel(siz) == 2 && all(siz == fix(siz) & siz >= 1))
    error('mrf2d:invalidarg', ...
        'The first argument [m, n] should be a pair of positive integer scalars.');
end

h = siz(1);
w = siz(2);

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

ss = cell(1, nk);
ts = cell(1, nk);
ws = cell(1, nk);

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
        [s, t] = make_links(w, h, x0, x0, y0, y1);
        
        ss{k} = s;
        ts{k} = t;
        
    elseif dy == 0
        [s, t] = make_links(w, h, x0, x1, y0, y0);
        
        ss{k} = s;
        ts{k} = t;
        
    else
        [s0, t0] = make_links(w, h, x0, x1, y0, y1);
        [s1, t1] = make_links(w, h, x0, x1, y1, y0);
        
        ss{k} = [s0; s1];
        ts{k} = [t0; t1];
    end

    ws{k} = constmat(numel(ss{k}), 1, kv);
end
    
s = vertcat(ss{:});
t = vertcat(ts{:});
w = vertcat(ws{:});

 

%% subfunctions

function [i, j] = make_links(w, h, xi, xj, yi, yj) %#ok<INUSL>

i = yi(:) + h * (xi(:) - 1);
j = yj(:) + h * (xj(:) - 1);


