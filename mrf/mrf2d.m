function W = mrf2d(siz, kernel, roi)
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
%   W = mrf2d([m, n], kernel, roi);
%       returns the affinity matrix of MRF model between the nodes
%       marked by roi. The nodes are re-indexed following the order
%       of find(roi).
%

% Created by Dahua Lin, on Apr 16, 2010
% Modified by Dahua Lin, on Apr 17, 2010
%   - add ROI support

%% verify input arguments

if ~(isnumeric(siz) && numel(siz) == 2 && all(siz == fix(siz) & siz >= 1))
    error('mrf2d:invalidarg', ...
        'The first argument [m, n] should be a pair of positive integer scalars.');
end

m = siz(1);
n = siz(2);

if ~(isfloat(kernel) && isreal(kernel) && ndims(kernel) == 2)
    error('mrf2d:invalidarg', 'kernel should be a real matrix.');
end

if kernel(1) ~= 0
    error('mrf2d:invalidarg', 'kernel(1, 1) must be zero.');
end

h = size(kernel, 1) - 1;
w = size(kernel, 2) - 1;

if m <= h || n <= w
    error('mrf2d:invalidarg', 'm and n should be greater than the kernel size.');
end

if nargin < 4
    roi = [];
    N = m * n;  % total number of nodes
else
    if ~(islogical(roi) && isequal(size(roi), [m n]))
        error('mrf2d:invalidarg', 'roi should be a logical matrix of size m x n.');
    end
    N = nnz(roi);
    imap = zeros(m, n);
    imap(roi) = 1 : N;
end


%% main

if h == 0 && w == 0  % no inter-node links
    
    W = sparse([], [], [], N, N);
        
else
    
    Is = cell(h+1, w+1);
    Js = cell(h+1, w+1);
    Vs = cell(h+1, w+1);
    
    for k = 0 : h
        for l = 0 : w
            kv = kernel(k+1, l+1);
            
            if k == 0 && l == 0
                continue;
                
            elseif k == 0
                [I1, J1, V1] = make_links(m, n, l, 0, kv);
                I = [I1; J1];
                J = [J1; I1];
                V = [V1; V1];
                    
            elseif l == 0
                [I1, J1, V1] = make_links(m, n, 0, k, kv);
                I = [I1; J1];
                J = [J1; I1];
                V = [V1; V1];
                
            else
                [I1, J1, V1] = make_links(m, n, l, k, kv);
                [I2, J2, V2] = make_links(m, n, l, -k, kv);
                I = [I1; J1; I2; J2];
                J = [J1; I1; J2; I2];
                V = [V1; V1; V2; V2];            
            
            end
            
            Is{k+1, l+1} = I;
            Js{k+1, l+1} = J;
            Vs{k+1, l+1} = V;
        end
    end
    
    I = vertcat(Is{:});
    J = vertcat(Js{:});
    V = vertcat(Vs{:});
    
    if isempty(roi)
        W = sparse(I, J, V, N, N);
    else
        si = roi(I) & roi(J);
        I = imap(I(si));
        J = imap(J(si));
        V = V(si);
        W = sparse(I, J, V, N, N);
    end
        
end
    

%% sub functions

function [I, J, V] = make_links(m, n, dx, dy, v)

[x0, x1] = shift_rows(n, dx);
[y0, y1] = shift_rows(m, dy);

[X0, Y0] = meshgrid(x0, y0);
[X1, Y1] = meshgrid(x1, y1);

I = Y0(:) + (X0(:) - 1) * m;
J = Y1(:) + (X1(:) - 1) * m;
V = v(ones(size(I)));



function [x0, x1] = shift_rows(n, d)

if d == 0
    x0 = 1 : n;
    x1 = x0;
elseif d > 0
    x0 = 1 : n - d;
    x1 = 1 + d : n;
else
    x0 = 1 - d : n;
    x1 = 1 : n + d;
end

    






