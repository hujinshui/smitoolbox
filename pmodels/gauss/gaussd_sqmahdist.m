function D = gaussd_sqmahdist(G, X, ca)
% Evaluate the squared Mahalanobis distances to Gaussian centers
%
%   D = gaussd_sqmahdist(G, X);
%   D = gaussd_sqmahdist(G, X, ca)
%
%       computes the Gaussian Mahalanobis distance to the centers of 
%       Gaussian distribution in G.
%
%       Given a Gaussian distribution with mean mu and covariance C,
%       the squared Mahalanobis distance of x w.r.t. this model is
%
%               (x - mu)' * inv(C) * (x - mu)
%
%       Inputs:
%       - G:        a gaussd struct.
%       - X:        the matrix comprised of samples as columns: size d x n.
%       - ca:       the value of mu' * C * mu. If ca is input, the input
%                   value will be used, otherwise, the function will call
%                   gaussd_const to compute it.
%
%       Outputs: (suppose K = G.n, and n is #samples in X)
%       - R:        a K x n matrix, where R(k, i) is the squared
%                   Mahalanobis distance of X(:,i) to the k-th model in
%                   G.
%

% Created by Dahua Lin, on Dec 5, 2011
%

%% verify inputs

if ~is_gaussd(G)
    error('gaussd_sqmahdist:invalidarg', 'G should be a gaussd struct.'); 
end

d = G.d;
if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == d)
    error('gaussd_sqmahdist:invalidarg', ...
        'X should be a real matrix with size(X,1) == d.');
end

if nargin < 3
    ca = [];
end


%% main

% type-specific computation of individual terms

ty = G.ty;

if ty == 'm'
    
    C = G.C;
    mu = G.mu;
    
    if isequal(mu, 0)
        zm = 1;
    else
        zm = 0;
        h = pdmat_lsolve(C, mu);
        M1 = h' * X;
    end

    if d == 1
        M2 = (1 ./ C.v(:)) * (X.^2);
    else
        switch C.ty
            case 's'
                M2 = (1 ./ C.v)' * sum(X.^2, 1);
            case 'd'
                M2 = (1 ./ C.v)' * (X.^2);
            case 'f'
                Cmat = C.v;
                if C.n == 1
                    M2 = dot(X, Cmat \ X, 1);
                else
                    M2 = zeros(C.n, size(X,2));
                    for k = 1 : C.n
                        M2(k,:) = dot(X, Cmat(:,:,k) \ X, 1);
                    end
                end
        end
    end
 
    
elseif ty == 'c'
    
    J = G.J;
    h = G.h;
    
    if isequal(h, 0)
        zm = 1;
    else
        zm = 0;
        M1 = h' * X;
    end    
    
    if d == 1
        M2 = J.v(:) * (X.^2);
    else
        switch J.ty
            case 's'
                M2 = J.v' * sum(X.^2, 1);
            case 'd'
                M2 = J.v' * (X.^2);
            case 'f'
                Jmat = J.v;
                if J.n == 1
                    M2 = dot(X, Jmat * X, 1);
                else
                    M2 = zeros(J.n, size(X,2));
                    for k = 1 : J.n
                        M2(k,:) = dot(X, Jmat(:,:,k) * X, 1);
                    end
                end
        end
    end
    
end

% combine terms

if zm
    D = M2;
else
    K = G.n;
    if isempty(ca)
        ca = gaussd_const(G);
    else        
        if ~(isfloat(ca) && isreal(ca) && isvector(ca) && numel(ca) == K)
            error('gaussd_sqmahdist:invalidarg', ...
                'The input ca should be a vector of length K.');
        end
    end
    
    if size(M1, 1) == size(M2, 1)
        D = M2 - 2 * M1;
    else
        D = bsxfun(@minus, M2, 2 * M1);
    end
    
    D = bsxfun(@plus, D, ca(:));    
end

D(D < 0) = 0;


