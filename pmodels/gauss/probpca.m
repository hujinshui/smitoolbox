classdef probpca
    % The class to represent a Probabilistic PCA model
    %
    
    % Created by Dahua Lin, on Nov 20, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        dim;        % the dimension of the (observed) vector space (d)
        ldim;       % the dimension of the latent space (q)
        
        mu;         % the mean vector (a vector of size d x 1 or a scalar 0)
        basis;      % the basis of the principal subspace (d x q)
        eigvals;    % the principal eigenvalues (q x 1)
        sigma2;     % the noise variance (scalar)
        
        W;          % the transform from latent to observed (d x q)
        T;          % the transform from observed to latent (q x d)
        
        
        % the following quantities are for evaluation of logpdf
        
        F;          % F := (sigma^2 I + W'W)^{-1} * W' (q x d)
                    % The information matrix has J = (1/sigma2) * (I - WF)                    
        h;          % the potential vector (d x 1)
        c0;         % the constant part: mu' * J * mu = mu' * h
        ldcov;      % the log-determinant of covariance        
    end
    
    
    methods
    
        function obj = probpca(mu, U, eigvals, sigma2)
            % Constructs a probabilistic PCA model
            %
            %   obj = probpca(mu, U, eigvals, sigma2);
            %       constructs a probabilistic PCA model.
            %
            %       Input arguments:
            %       - mu:       the mean vector (d x 1 vector or just 0)
            %       - U:        the basis of principal subspace (d x q)
            %       - eigvals:  the principal eigenvalues (q x 1)
            %       - sigma2:   the noise variance (scalar)            
            %
            
            % verify input
            
            if ~(isfloat(U) && ndims(U) == 2)
                error('probpca:invalidarg', ...
                    'U should be a numeric matrix.');
            end
            [d, q] = size(U);
            
            if ~(isfloat(mu) && (isequal(size(mu), [d 1]) || isequal(mu, 0)))
                error('probpca:invalidarg', ...
                    'mu should be either a d x 1 vector or a scalar.');
            end
            
            if ~(isfloat(eigvals) && isvector(eigvals) && numel(eigvals) == q)
                error('probpca:invalidarg', ...
                    'eigvals should be a numeric vector of length q.');
            end
            if size(eigvals, 2) > 1; eigvals = eigvals.'; end            
            
            if ~(isfloat(eigvals) && isscalar(sigma2) && sigma2 > 0)
                error('probpca:invalidarg', ...
                    'sigma2 should be a positive scalar.');
            end
            
            
            % set basic fields
            
            obj.dim = d;
            obj.ldim = q;
            obj.mu = mu;
            obj.basis = U;
            obj.eigvals = eigvals;
            obj.sigma2 = sigma2;
            
            % compute W and T
            
            s = sqrt(eigvals - sigma2);
            W_ = bsxfun(@times, U, s.');
            T_ = bsxfun(@times, U', 1./s);
            
            obj.W = W_;
            obj.T = T_;
            
            % compute information parameters
                            
            F_ = adddiag(W_' * W_, sigma2) \ W_';
            if isequal(mu, 0)
                h_ = 0;
                c0_ = 0;
            else
                h_ = (mu - W_ * (F_ * mu)) * (1 / sigma2);
                c0_ = h_' * mu;
            end
            
            ldc = sum(log(eigvals)) + (d - q) * log(sigma2);
            
            obj.F = F_;
            obj.h = h_;
            obj.c0 = c0_;
            obj.ldcov = ldc;
            
        end
        
        
        function C = get_cov(obj)
            % Compute the covariance matrix
            %
            %   C = obj.get_cov();
            %       It returns a d x d covariance matrix.
            %
            %       Note that C is a dense matrix, and be cautious in
            %       using this method when d is very large.
            %
            
            W_ = obj.W;
            C = adddiag(W_ * W_', obj.sigma2);
        end   
        
        
        function G = to_gauss(obj, op)
            % Convert the PPCA model to a gaussd object
            %
            %   G = to_gauss(obj, 'mp');
            %   G = to_gauss(obj, 'ip');
            %   G = to_gauss(obj, 'both');
            %       It returns an mathematically equivalent gaussd
            %       object, with covariance represented by a gsymat
            %       object.
            %
            %       Note that this function involves explicitly
            %       computing the full covariance matrix or full 
            %       information matrix. Be cautious
            %       in using this function when d is very large.
            %
            %       The returned object have mean parameter and/or
            %       information parameters, depending on the 2nd arg.
            %
            
            if ~ischar(op) 
                error('probpca:to_gauss:invalidarg', ...
                    'The 2nd argument should be a string.');
            end
            
            switch op
                case 'mp'
                    C = get_cov(obj);
                    G = gaussd.from_mp(obj.mu, gsymat(C));
                case 'ip'
                    J = obj.W * obj.F;
                    J = (-0.5) * (J + J');
                    J = (1/obj.sigma2) * adddiag(J, 1);
                    G = gaussd.from_ip(obj.h, gsymat(J), obj.c0);
                case 'both'
                    C = get_cov(obj);
                    G = gaussd.from_mp(obj.mu, gsymat(C), 'ip');
                otherwise
                    error('probpca:invalidarg', ...
                        'The 2nd argument is invalid.');
            end                                
        end
        
        
        function X = reconstruct(obj, Y)
            % reconstructs from latent variables
            %
            %   X = obj.reconstruct(Y);
            %       reconstructs the vector in the d-dimensional space
            %       based on the latent representation in Y, as follows.
            %
            %           x = W y + mu
            %   
            %       In the output, X(:,i) corresponds to Y(:,i).
            %
            
            X = obj.W * Y;
            mu_ = obj.mu;
            if ~isequal(mu_, 0)
                X = bsxfun(@plus, X, mu_);
            end
        end
        
        
        function Y = transform(obj, X)
            % transform observed vectors to latent representation.
            %
            %   Y = obj.transform(X);
            %       transforms observed vectors to the q-dimensional
            %       representation on latent space.
            %
            
            mu_ = obj.mu;
            if ~isequal(mu_, 0)
                X = bsxfun(@minus, X, mu_);
            end            
            Y = obj.T * X;          
        end
       
        
        function D = sqmahdist(obj, X)
            % Compute the Mahalanobis distance of X to center
            %
            %   D = obj.sqmahdist(X);
            %
            %       X is the input sample matrix, whose size can be a 
            %       d x n matrix, and then D will be a row vector of 
            %       length n, with D(i) corresponding to X(:,i).
            %
            
            d = obj.dim;
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == d)
                error('probpca:sqmahdist:invalidarg', ...
                    'X should be a numeric matrix with size(X,1) == d.');
            end
            
            h_ = obj.h;
            W_ = obj.W;
            F_ = obj.F;
            
            Z = X - W_ * (F_ * X);
            t2 = dot(X, Z, 1) * (1 / obj.sigma2);
            
            if isequal(h_, 0)
                D = t2;
            else
                c0_ = obj.c0;
                t1 = h_' * X;
                
                D = t2 - 2 * t1 + c0_;
            end                
        end
        
        
        function L = logpdf(obj, X)
            % Compute the log-pdf of observed vectors
            %
            %   L = obj.logpdf(X);
            %       Suppose there are n columns in X, then L will be
            %       a row vector of length n. In particular, L(i)
            %       corresponds to X(:,i).
            %
            
            D = sqmahdist(obj, X);
            L = -0.5 * (D + (obj.dim * log(2*pi) + obj.ldcov));
        end
                
    end
    
    
    methods(Static)
       
        function obj = from(varargin)
            % Constructs a PPCA model from covariance or standard PCA
            %
            %   obj = probpca.from(mu, C, p);
            %       constructs a Probabilistic PCA model with latent 
            %       dimension p from a Gaussian distribution, whose
            %       mean is mu, and covariance is C.
            %
            %       Note that p should be less than rank(C).
            %
            %   obj = probpca.from(pca);
            %   obj = probpca.from(pca, p);
            %       constructs a PPCA model with latetn dimension p
            %       from a standard PCA model. Here, pca should be
            %       an object of class pca_std.
            %
            %       When p is omitted, it uses pca.pdim as p.
            %
            
            if isnumeric(varargin{1}) % from mu and cov
                
                % verify argument
                
                u = varargin{1};
                C = varargin{2};
                p = varargin{3};
                                
                d = size(C, 1);
                if ~(isfloat(C) && ndims(C) == 2 && d == size(C,2))
                    error('probpca:from:invalidarg', ...
                        'C should be a d x d matrix.');
                end
                
                if ~(isfloat(u) && (isequal(u, 0) || isequal(size(u), [d 1])))
                    error('probpca:from:invalidarg', ...
                        'mu should be either 0 or a column vector of size d x 1.');
                end                
                
                if ~(isnumeric(p) && isscalar(p) && p == fix(p) && p >= 1 && p < d)
                    error('probpca:from:invalidarg', ...
                        'p should be an integer in [1, d-1]');
                end
                
                % compute
                
                [U, evs] = eig(C);
                evs = diag(evs);
                
                [evs, si] = sort(evs, 1, 'descend');
                pevs = evs(1:p);
                Up = U(:, si(1:p));
                
                sig2 = sum(evs(p+1:end)) / (d-p);
                
                % construct 
                
                obj = probpca(u, Up, pevs, sig2);                
                
            elseif isa(varargin{1}, 'pca_std')
                
                % verify argument
                
                S = varargin{1};
                sp = S.pdim;
                
                if nargin < 2
                    p = sp;
                else
                    p = varargin{2};
                    if ~(isnumeric(p) && isscalar(p) && p == fix(p) && p >= 1 && p <= sp)
                        error('probpca:from:invalidarg', ...
                            'p should be an integer in [1, pca.dim]');
                    end
                end
                
                % make
                
                if p == sp
                    Up = S.basis;
                    pevs = S.eigvals;
                    sig2 = S.residue_var / (S.dim - p);
                else
                    Up = S.basis(:, 1:p);
                    pevs = S.eigvals(1:p);
                    sig2 = (sum(S.eigvals(p+1:end)) + S.residue_var) / (S.dim - p);
                end
                
                obj = probpca(S.center, Up, pevs, sig2);                
            end            
        end
        
                
    end
    
        
end


