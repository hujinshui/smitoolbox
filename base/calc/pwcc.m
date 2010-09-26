function V = pwcc(f, X, Y)
% Perform pairwise computation between columns
%
%   V = pwcc(f, X, Y);
%       applies the function f to columns in X and Y pairwisely.
%
%       f should be a function that can support the following syntax:
%       
%           v = f(A, B);
%
%       where A and B have the same number of columns n, then v has size
%       1 x n, such that v(i) corresponds to A(:,i) and B(:,i).
%       Note that the number of rows in A and B need not be equal. 
%
%       In the output, V has size n1 x n2, where n1 and n2 are respectively
%       the number of columns in X and Y.

%   Created by Dahua Lin, on Sep 12, 2010
%

%% verify input

if ndims(X) > 2 || ndims(Y) > 2
    error('pwcc:invalidarg', 'X and Y should be matrices.');
end

%% main

nx = size(X, 2);
ny = size(Y, 2);

if nx <= ny    
    if nx == 1        
        V = f(X(:, ones(1, ny)), Y);        
    else
        x = X(:, 1);        
        v = f(x(:, ones(1, ny)), Y);
        V = zeros(nx, ny, class(v));
        V(1, :) = v;
        for i = 2 : nx
            x = X(:, i);
            v = f(x(:, ones(1, ny)), Y);
            V(i, :) = v;
        end
    end    
else     
    if ny == 1               
        V = f(X, Y(:, ones(1, nx))).';        
    else
        y = Y(:, 1);
        v = f(X, y(:, ones(1, nx)));
        v = zeros(nx, ny, class(v));
        V(:, 1) = v.';
        for i = 2 : ny
            y = Y(:, i);
            v = f(X, y(:, ones(1, nx)));
            V(:, i) = v.';
        end
    end    
end
    
