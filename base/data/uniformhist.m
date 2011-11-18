function [H, varargout] = uniformhist(X, hsiz, rgn)
% Generate histogram over uniform bins
%
%   [H, C] = uniformhist(X, n, [a, b]);
%       
%       generates a histogram over uniform bins for one-dimensional
%       data. Here, X is a matrix comprised of observed values. 
%       
%       The interval [a, b] is divided into n uniform bins.
%       
%       In the output, H is a vector of size n x 1, in which H(i) is the
%       number of values in the i-th bin. 
%
%       C is also a vector of n x 1, where C(i) is the center
%       of the i-th bin.
%
%   [H, C1, C2] = uniformdist(X, [m, n], [a1, b1, a2, b2]);
%
%       generates a histogram over uniform bins for two-dimensional
%       points. Here, X should be a matrix of size N x 2, with
%       each row being a point.
%
%       The interval [a1, b1] is divided into m uniform bins, and
%       the interval [a2, b2] is divided into n uniform bins. 
%       In this way, we obain m x n rectangle bins over a 2D space.
%
%       In the output, H is a matrix of size m x n, in which H(i,j) 
%       is the number of points in the bin at i-th row and j-th column.
%
%       centers is a cell arrray, where centers{1} is a vector of size
%       m x 1, where C1(i) is the center of the i-th bin along the first 
%       dimension, and C2(j) is the center of the j-th bin along the 
%       second dimension.
%
%   [H, C1, ..., Cd] = uniformdist(X, hsiz, rgn);
%
%       This function can be applied to the data of arbitrary dimensions.
%       
%       Input arguments:
%       - X:        an N x d data matrix, each row is a point
%       - rgn:      the range, in form of [a1, b1, ..., ad, bd]
%       - hsiz:     the size of histogram, in form of [n1, ..., nd]
%

% Created by Dahua Lin, on Nov 2, 2011
%

%% veirfy input arguments

if ~(isfloat(X) && ndims(X) == 2 && isreal(X) && ~issparse(X))
    error('uniformhist:invalidarg', 'X should be a non-sparse real matrix.');
end

if ~(isnumeric(hsiz) && ndims(hsiz) == 2 && size(hsiz, 1) == 1)
    error('uniformhist:invalidarg', 'hsiz should be a numeric row vector.');
end
d = numel(hsiz);

if ~(isnumeric(rgn) && isequal(size(rgn), [1 2*d]))
    error('uniformhist:invalidarg', 'rgn is invalid.');
end

if d == 1
    if ~isvector(X)
        X = X(:);
    end
else
    if size(X,2) ~= d
        error('uniformhist:invalidarg', 'X should have %d columns.', d);
    end
end
    

%% main

if d == 1
    
    a = rgn(1);
    b = rgn(2);
    n = hsiz;
    
    s = (b - a) / n;
    Z = ceil((1 / s) * (X - a));
    H = intcount(n, Z)';
        
    C = a + (0.5 + (0:n-1)) * s;  
    varargout{1} = C;
    
elseif d == 2
    
    a1 = rgn(1);
    b1 = rgn(2);
    a2 = rgn(3);
    b2 = rgn(4);
    
    n1 = hsiz(1);
    n2 = hsiz(2);
    
    s1 = (b1 - a1) / n1;
    s2 = (b2 - a2) / n2;
    
    x1 = X(:,1);
    x2 = X(:,2);
    
    r = find(x1 >= a2 & x1 <= b2 & x2 >= a2 & x2 <= b2);
    x1 = x1(r);
    x2 = x2(r);
    
    Z1 = ceil((1 / s1) * (x1 - a1));
    Z2 = ceil((1 / s2) * (x2 - a2));
    
    Z = Z1 + (Z2 - 1) * n1;
    H = intcount(n1 * n2, Z);
    H = reshape(H, n1, n2);
    
    C1 = a1 + (0.5 + (0:n1-1)) * s1;
    C2 = a2 + (0.5 + (0:n2-1)) * s2;
    
    varargout = {C1, C2};

else % d > 2
    
    a = rgn(1:2:end);
    b = rgn(2:2:end);
    
    N = size(X, 1);
    
    % filter data
    
    r = true(N, 1);
    for k = 1 : d
        x = X(:, k);
        r = r & (x >= a(k) & x <= b(k));
    end
    
    X = X(r, :);
    
    % compute H
    
    Zs = cell(1, d);
    for k = 1 : d
        x = X(:, k);
        s = (b(k) - a(k)) / hsiz(k);
        z = ceil((1 / s) * (x - a(k)));
        Zs{k} = z;
    end
    
    Z = sub2ind(hsiz, Zs{:});
    H = intcount(prod(hsiz), Z);   
    H = reshape(H, hsiz);
    
    % compute C
    varargout = cell(1, d);
    for k = 1 : d
        cn = hsiz(k);
        s = (b(k) - a(k)) / cn;
        c = a(k) + (0.5 + (0:cn-1)) * s;
        varargout{k} = c;
    end
    
end


