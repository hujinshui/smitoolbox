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

% compute M2 (and at the same time mu or h)

ty = G.ty;

if ty == 'm'    % with mean params
    
    C = G.C;
    mu = G.mu;
    
    if isequal(mu, 0)
        zm = 1;
    else
        zm = 0;
    end
        
    if d == 1
        r = 1 ./ C.v(:);
        M2 = r * (X.^2);
        
        if ~zm
            if G.n == 1
                h = r * mu;
            else
                h = bsxfun(@times, r.', mu);
            end
        end        
    else
        switch C.ty
            case 's'
                r = 1 ./ C.v;
                M2 = r' * sum(X.^2, 1);
                
                if ~zm
                    if G.n == 1
                        h = r * mu;
                    else
                        h = bsxfun(@times, r, mu);
                    end
                end
                
            case 'd'
                r = 1 ./ C.v;
                M2 = r' * (X.^2);
                
                if ~zm
                    if C.n == G.n
                        h = r .* mu;
                    else
                        h = bsxfun(@times, r, mu);
                    end
                end
                
            case 'f'
                Cm = C.v;
                if C.n == 1                    
                    if zm
                        JX = Cm \ X;
                    else
                        [JX, h] = solve_two(Cm, X, mu);
                    end
                    
                    M2 = dot(X, JX, 1);

                else
                    M2 = zeros(C.n, size(X,2));
                    
                    if zm                    
                        for k = 1 : C.n
                            M2(k, :) = dot(X, Cm(:,:,k) \ X, 1);
                        end
                    else
                        h = zeros(d, G.n);
                        for k = 1 : C.n
                            [JX, h(:,k)] = solve_two(Cm(:,:,k), X, mu(:,k));
                            M2(k, :) = dot(X, JX, 1);
                        end
                    end
                end
        end
    end    
    
elseif ty == 'c'    % with canonical params
    
    J = G.J;
    h = G.h;
    
    if isequal(h, 0)
        zm = 1;
    else
        zm = 0;
    end    
    
    M2 = pdmat_quad(J, X);
    
    if isempty(ca) && ~zm      % need mu to calculate ca
        mu = pdmat_lsolve(J, h);
    end
end

% combine terms

if zm
    D = M2;
else
    K = G.n;
    if isempty(ca)
        if K == 1
            ca = h' * mu;
        else
            ca = dot(h, mu, 1);            
        end
    else        
        if ~(isfloat(ca) && isreal(ca) && isvector(ca) && numel(ca) == K)
            error('gaussd_sqmahdist:invalidarg', ...
                'The input ca should be a vector of length K.');
        end
    end
    
    M1 = h' * X;
    if size(M1, 1) == size(M2, 1)
        D = M2 - 2 * M1;
    else
        D = bsxfun(@minus, M2, 2 * M1);
    end
    
    D = bsxfun(@plus, D, ca(:));    
end

D(D < 0) = 0;


%% auxiliary function

function [JX, h] = solve_two(C, X, mu)

n = size(X, 2);
JX = C \ [X mu];
h = JX(:, n+1:end);
JX = JX(:, 1:n);

