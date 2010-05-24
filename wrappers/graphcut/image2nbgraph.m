function G = image2nbgraph(I, wfun, kernel, roi)
%IMAGE2NBGRAPH Constructs a neighborhood graph from an image
%
%   G = image2nbgraph(I, wfun, kernel, roi);
%
%   Input:
%       I:      the input image [a 2D or 3D numeric array]
%       wfun:   the function to compute weights between two lists of pixels
%               wfun can be empty, then it assigns weight 1 for every pair.
%       kernel: the coefficient kernel to modulate the weights
%
%   Output:
%       G:      the constructed graph in form of a struct with following fields:
%               - n:        the number of nodes
%               - m:        the number of edges
%               - edges:    the edges [2 x m matrix]
%               - eweights: the weights of the edges [1 x m vector]
%

%% parse and verify input arguments

error(nargchk(3, 3, nargin));

assert(isnumeric(I) && ndims(I) <= 3, 'imag2nbgraph:invalidarg', ...
    'The image I should be a 2D or 3D numeric array.');

h = size(I, 1);
w = size(I, 2);

if ~isempty(wfun)
    assert(isa(wfun, 'function_handle'), 'image2nbgraph:invalidarg', ...
        'wfun should be a function handle.');
end

assert(isfloat(kernel) && ndims(kernel) == 2, 'image2nbgraph:invalidarg', ...
    'kernel should be a numeric matrix.');

%% main


[dys, dxs, kvs] = find(kernel);
dxs = dxs - 1;
dys = dys - 1;
nkv = numel(kvs);

if ~isempty(wfun)
    if ndims(I) == 2
        V = reshape(I, 1, h * w);
    else
        pd = size(I, 3);
        V = reshape(permute(I, [3, 1, 2]), [pd, h * w]);
    end
end
    

esets = cell(nkv, 1);
ewsets = cell(nkv, 1);

for i = 1 : nkv
    
    dx = dxs(i);
    dy = dys(i);    
    kv = kvs(i);   
    
    % extract neighboring links
    
    nlinks = get_nlinks(w, h, dx, dy);
    
    if isempty(nlinks)
        continue;
    end    
    
    % compute edge weights
    
    if ~isempty(wfun)
        ews = wfun(V(:, nlinks(1, :)), V(:, nlinks(2, :))) * kv;
            
        assert(isnumeric(ews) && isequal(size(ews), [1, size(nlinks, 2)]), ...
            'image2nbgraph:rterror', 'the result by wfun is invalid.');
    else
        ews = ones(1, size(nlinks, 2));
    end
        
    
    % store
    
    esets{i} = nlinks;
    ewsets{i} = ews;    
end

% output

edges = [esets{:}];
eweights = [ewsets{:}];

G = struct( ...
    'n', h * w, ...
    'm', size(edges, 2), ...
    'edges', edges, ...
    'eweights', eweights);


function nlinks = get_nlinks(w, h, dx, dy)
% The function to extract neighbor links between pairs of pixels of specified 
% displacement

nlinks = [];

if w <= dx || h <= dy || (dx == 0 && dy == 0)
    return;
end

[x1, y1] = meshgrid(1 : w-dx, 1 : h-dy);

x1 = reshape(x1, 1, []);
y1 = reshape(y1, 1, []);

x2 = x1 + dx;
y2 = y1 + dy;

idx1 = sub2ind([h, w], y1, x1);
idx2 = sub2ind([h, w], y2, x2);

nlinks = [idx1; idx2];




