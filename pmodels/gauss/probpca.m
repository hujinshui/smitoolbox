classdef probpca
    % The class to represent a Probabilistic PCA model
    %
    %   A Probabilistic PCA (PPCA) model is parameterized by a d x q
    %   matrix W that relates the latent variables with the observed
    %   variables, and sigma^2.
    %
    %       z ~ N(0, I);                    -- the latent variable x
    %       x ~ N(W * z + mu, sigma^2 * I)  -- the observed vector y
    %
    
    % Created by Dahua Lin, on Nov 20, 2010
    % Modified by Dahua Lin, on Nov 3, 2011
    %
    
    %% Properties
    
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
    
    
    %% Constructor
    
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
    end
    
    
    %% Basic methods
    
    methods        
        
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
        
        
        function J = get_precmat(obj)
            % Compute the precision matrix (inverse covariance)
            %
            %   J = obj.get_precmat();
            %       It returns a d x d precision matrix.
            %
            %       Note that C is a dense matrix, and be cautious in
            %       using this method when d is very large.
            %
            
            J = obj.W * obj.F;
            J = (-0.5) * (J + J');
            J = (1/obj.sigma2) * adddiag(J, 1);            
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
                    G = gaussd.from_mp(obj.mu, pdmat(C));
                case 'ip'
                    J = get_precmat(obj);
                    G = gaussd.from_ip(obj.h, pdmat(J), obj.c0);
                case 'both'
                    C = get_cov(obj);
                    G = gaussd.from_mp(obj.mu, pdmat(C), 'ip');
                otherwise
                    error('probpca:invalidarg', ...
                        'The 2nd argument is invalid.');
            end                                
        end
        
        
        function X = transform_fwd(obj, Z)
            % reconstructs from latent variables
            %
            %   X = obj.transform_fwd(Y);
            %       reconstructs the vector in the d-dimensional space
            %       based on the latent representation in Y, as follows.
            %
            %           x = W y + mu
            %   
            %       In the output, X(:,i) corresponds to Y(:,i).
            %
            
            X = obj.W * Z;
            mu_ = obj.mu;
            if ~isequal(mu_, 0)
                X = bsxfun(@plus, X, mu_);
            end
        end
        
        
        function Z = transform_bwd(obj, X)
            % transform observed vectors to latent representation.
            %
            %   Z = obj.transform_bwd(X);
            %       transforms observed vectors to the q-dimensional
            %       representation on latent space.
            %
            
            mu_ = obj.mu;
            if ~isequal(mu_, 0)
                X = bsxfun(@minus, X, mu_);
            end            
            Z = obj.T * X;          
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
            t2 = sum(X .* Z, 1) * (1 / obj.sigma2);
            
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
            
        
        function P = pdf(obj, X)
            % Compute the pdf of observed vectors in X
            %
            %   L = obj.pdf(X);
            %       Suppose there are n columns in X, then P will be
            %       a row vector of length n. In particular, P(i)
            %       corresponds to X(:,i).
            %
            
            P = exp(logpdf(obj, X));
        end
        
        
        function X = sample(obj, n)
            % Draw samples from the probabilistic PCA model
            %
            %   X = sample(obj, n);
            %
            %       Draws n samples from the PPCA model. 
            %       In the output, X will be a d x n matrix.
            %
            
            d = obj.dim;
            q = obj.ldim;
            
            Z = randn(q, n);
            X = obj.W * Z;
            
            u = obj.mu;
            if ~isequal(u, 0)
                X = bsxfun(@plus, X, u);
            end
            
            s = sqrt(obj.sigma2);
            X = X + randn(d, n) * s;            
        end
        
    end
    
    
    %% Estimation
    
    methods(Static)
        
        function obj = mle(X, w, q, varargin)
            % Performs Maximum Likelihood estimation of PPCA from data
            %
            %   obj = probpca.mle(X, [], q, ...);
            %   obj = probpca.mle(X, w, q, ...);
            %
            %       performs Maximum Likelihood estimation of the PPCA
            %       model from data. 
            %
            %       Input arguments:
            %       - X:    The data matrix of size d x n, of which each
            %               column is s sample.
            %       - w:    The sample weights, a row vector of size 
            %               1 x n. If all samples have the same weight,
            %               it can be empty.
            %       - q:    the dimension of the latent space. It should 
            %               have q < min(d, n).
            %
            %       One can specify further options to control the
            %       estimation, in form of name/value pairs.
            %
            %       - 'method':     The method used to do the training,
            %                       which can be
            %                       - 'cov':    by computing the covariance
            %                                   matrix, and compute the
            %                                   eigenvectors of it.
            %                       - 'std':    by doing SVD, this can be
            %                                   faster when n < d.
            %                       The default is 'cov'.
            %
            %       - 'zmean':      Whether to assume a zero mean vector.
            %
            
            % verify arguments
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
                error('probpca:mle:invalidarg', ...
                    'X should be a real matrix.');
            end
            [d, n] = size(X);
            
            if ~isempty(w)
                if ~(isfloat(w) && isequal(size(w), [1 n]))
                    error('probpca:mle:invalidarg', ...
                        'w should be a vector of size 1 x n.');
                end
            end
            
            if ~(isnumeric(q) && isscalar(q) && q == fix(q) && ...
                    q >= 1 && q < d && q < n);
                error('probpca:mle:invalidarg', ...
                    'q should be a positive integer less than min(d, n).');
            end
            
            % check options
            
            method = 'cov';
            zmean = false;
            
            if ~isempty(varargin)
                onames = varargin(1:2:end);
                ovals = varargin(2:2:end);
                
                if ~(numel(onames) == numel(ovals) && iscellstr(onames))
                    error('probpca:mle:invalidarg', ...
                        'The name/value list for options is invalid.');
                end
                
                for i = 1 : numel(onames)
                    cn = onames{i};
                    cv = ovals{i};
                    switch lower(cn)
                        case 'method'
                            if ~(ischar(cv) && ...
                                    (strcmp(cv, 'cov') || strcmp(cv, 'std')))
                                error('probpca:mle:invalidarg', ...
                                    'The method should be either ''cov'' or ''std''.');
                            end
                            method = cv;
                        case 'zmean'
                            if ~(islogical(cv) && isscalar(cv))
                                error('probpca:mle:invalidarg', ...
                                    'zmean should be a logical scalar');
                            end
                            zmean = cv;
                        otherwise
                            error('probpca:mle:invalidarg', ...
                                'Unknown option name %s', cn);
                    end
                end
            end
            
            % do estimation
            
            if zmean
                u = 0;
                Z = u;
            else
                if isempty(w)
                    u = sum(X, 2) * (1 / n);
                else
                    sw = sum(w);
                    u = (X * w') * (1 / sw);
                end
                Z = bsxfun(@minus, X, u);
            end
            
            switch method
                case 'cov'
                    if isempty(w)
                        C = (Z * Z') * (1/n);
                    else
                        C = (Z * bsxfun(@times, Z, w)') * (1/sw);
                        C = 0.5 * (C + C');
                    end
                    [U, evs] = eig(C);
                    evs = diag(evs);                    
                        
                case 'svd'
                    if isempty(w)
                        [U, svs] = svd(Z, 0); 
                        svs = diag(svs);
                        evs = (svs.^2) * (1/n);
                    else
                        [U, svs] = svd(bsxfun(@times, Z, w), 0);
                        svs = diag(svs);
                        evs = (svs.^2) * (1/sw);
                    end                                        
            end
            
            obj = probpca.make_obj(u, U, evs, q);
        end
        
        
        function obj = from_stats(mu, C, q)
            % Constructs a PPCA model from a statistics
            %
            %   obj = probpca.from(mu, C, q);
            %       constructs a PPCA model with latent dimension q, 
            %       from a mean vector mu and a covariance matrix C.
            %
            %       Suppose the space of observed vectors has dimension d.
            %       Then mu can be a d x 1 vector, and simply zero.
            %       C is a covariance matrix of size d x d.
            %
            
            % verify arguments
            
            if ~(isfloat(C) && ndims(C) == 2 && d == size(C,2))
                error('probpca:from_stats:invalidarg', ...
                    'C should be a d x d matrix.');
            end
            
            if ~(isfloat(mu) && (isequal(mu, 0) || isequal(size(mu), [d 1])))
                error('probpca:from_stats:invalidarg', ...
                    'mu should be either 0 or a column vector of size d x 1.');
            end
            
            if ~(isnumeric(q) && isscalar(q) && q == fix(q) && q >= 1 && q < d)
                error('probpca:from_stats:invalidarg', ...
                    'q should be an integer in [1, d-1]');
            end
            
            % compute
            
            [U, evs] = eig(C);
            evs = diag(evs);
            
            obj = probpca.make_obj(mu, U, evs, q);                        
        end                      
                
    end
        
    
    methods(Static, Access='private')
        
        function obj = make_obj(u, U, evs, q)
            
            d = size(U, 1);
            [evs, si] = sort(evs, 1, 'descend');
            qevs = evs(1:q);
            Uq = U(:, si(1:q));
            
            sig2 = sum(evs(q+1:end)) / (d-q);
            
            % construct
            
            obj = probpca(u, Uq, qevs, sig2);
            
        end        
        
    end
    
        
end


