function [s, t] = localedges(siz, rgn)
% Get the edges of 1D or 2D locally connected graph
%
%   [s, t] = localedges(n, r);
%
%       Gets the edges for locally connected graph over one-dimensional
%       grid. Here, n is the number of nodes, and r is the local range.
%       In this graph, each node is connected to a neighboring node,
%       within distance r.
%
%       In the output, s and t are vectors of size #E x 1, where #E is 
%       the number of undirected edges. In each pair of s(i) and t(i),
%       we have s(i) < t(i).
%   
%   [s, t] = localedges([m, n], r);
%   [s, t] = localedges([m, n], [yr, xr]);
%
%       Gets the edges for locally connected graph over two-dimensional
%       grid. Here, the grid size is m rows and n columns, and yr and xr
%       are respectively the local distance along y and x direction.
%
%       In the output, s and t are vectors of size #E x 2, where #E is 
%       the number of undirected edges. s(i,:) and t(i,:) are the
%       subscripts of the source and target end of the i-th edge
%
%   [s, t] = localedges([m, n], ker);
%
%       Gets the edges for locally connected graph over two-dimensional
%       grid. Here, ker is a matrix specifying the connection pattern, 
%       the node (i, j) is connected to (i+di, j+dj) if ker(di+1, dj+1) 
%       is true. 
%

%   Created by Dahua Lin, on Oct 27, 2011
%

%% verify input arguments

if ~(isnumeric(siz) && ndims(siz) == 2 && all(siz == fix(siz)) )
    error('localedges:invalidarg', 'The 1st argument is invalid.');
end

d = numel(siz);
if ~(d == 1 || d == 2)
    error('localedges:invalidarg', 'Only 1D or 2D grid is allowed.');
end

if d == 1
    r = rgn;
    if ~(isnumeric(r) && isscalar(r) && r == fix(r) && r >= 0)
        error('localedges:invalidarg', 'r should be a positive integer.');
    end
    
else
    if ndims(rgn) ~= 2
        error('localedges:invalidarg', 'The 2nd argument is invalid.');
    end
    
    if isnumeric(rgn) && numel(rgn) <= 2
        if isscalar(rgn)
            pat = true(rgn+1, rgn+1);
        else
            pat = true(rgn(1)+1, rgn(2)+1);
        end
        pat(1) = false;
        
    elseif islogical(rgn)        
        pat = rgn;         
        if pat(1)
            error('localedges:invalidarg', ...
                'The first element of the kernel must be false for 2D.');
        end
        
    else
        error('localedges:invalidarg', 'The 2nd argument is invalid.');
    end
    
end
    

%% main

if d == 1    
    n = siz;
    
    s = 1:n;
    s = s(ones(r, 1), :);
    t = bsxfun(@plus, (1:r).', 1:n);
    
    msk = t <= n;
    s = s(msk);
    t = t(msk);
        
else % d == 2
    
    h = siz(1);
    w = siz(2);
    
    % preprocess connection pattern
    
    kh = size(pat, 1);
    pat = pat([kh:-1:2, 1:kh], :);
    pat(1:kh) = false;
    
    [di, dj] = find(pat);
    di = di - kh;
    dj = dj - 1;
    
    ks = numel(di);
    
    % compute
    
    I = (1:h).';
    I = I(:, ones(1, w));
    I = reshape(I, 1, numel(I));
    
    J = 1:w;
    J = J(ones(h,1), :);
    J = reshape(J, 1, numel(J));
    
    I2 = bsxfun(@plus, I, di);
    J2 = bsxfun(@plus, J, dj);
    
    I = I(ones(ks, 1), :);
    J = J(ones(ks, 1), :);
    
    sel = find(I2 >= 1 & I2 <= h & J2 <= w);
    I = I(sel);
    J = J(sel);
    I2 = I2(sel);
    J2 = J2(sel);
    
    % combine
    
    s = [I J];
    t = [I2, J2];
        
end




