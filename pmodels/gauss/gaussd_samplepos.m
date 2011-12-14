function X = gaussd_samplepos(G, dh, dJ, n)
% Performs MAP estimation w.r.t. Gaussian prior
%
%   X = gaussd_samplepos(G, dh, dJ);
%   X = gaussd_samplepos(G, dh, dJ, n);
%
%       Samples from the posterior Gaussian distribution derived by
%       incorporating the given prior with the conjugate updates.
%       
%       Input arguments:
%       - G:    the Gaussian prior (G.ty == 'c' && G.n == 1)
%       - dh:   the update to the potential vector
%       - dJ:   the update to the precision matrix.
%
%       - n:    the number of samples to draw. If omitted, n is set to 1.
%               Note: n must be one, when size(dh, 2) > 1, meaning that
%               we are to draw from multiple posterior distribution. 
%               In this case, we draw one from each.
%
%       Output arguments:
%       - X:    the obtained samples.
%       

% Created by Dahua Lin, on Dec 14, 2011
%

%% verify input numbers

if nargin < 4
    n = 1;
end

m = size(dh, 2);
if m > 1 && n > 1
    error('gaussd_samplepos:invalidarg', ...
        'n must be one when size(dh, 2) > 1.');
end    

%% main

[h, J] = gaussd_conjupdate(G, dh, dJ);
d = size(h, 1);

if m == 1    
    X = gsample_c(randn(d, n), h, J);
else
    Z = randn(d, m);
    X = zeros(d, m);
        
    if J.n == 1        
        for i = 1 : m
            X(:,i) = gsample_c(Z(:,i), h(:,i), J);
        end
    else
        for i = 1 : m
            Ji = pdmat_pick(J, i);
            X(:,i) = gsample_c(Z(:,i), h(:,i), Ji);
        end
    end
end


%% core functions


function X = gsample_c(X, h, J)

[d, n] = size(X);
ty = J.ty;
v = J.v;

if ty == 's' || ty == 'd'
    
    if ~isequal(v, 1);
        v = 1 ./ v;
        if isequal(h, 0)
            mu = 0;
        else
            mu = h .* v;
        end

        if isscalar(v) || n == 1
            X = X .* sqrt(v);
        else
            X = bsxfun(@times, X, sqrt(v));
        end
    end
        
elseif ty == 'f'
    L = chol(v);
    
    if isequal(h, 0)
        mu = 0;
        X = L \ X;
    else                
        g = L' \ h;
        A = L \ [X g];
        X = A(:, 1:n);
        mu = A(:, n+1);
    end    
end

if ~isequal(mu, 0)
    if d == 1 || n == 1
        X = mu + X;
    else
        X = bsxfun(@plus, X, mu);
    end
end

