classdef gaussd
    % The class to represent Gaussian distribution(s)
    %
    %   Each object of this class can contain one or multiple Gaussian
    %   distributions. A Gaussian distribution can be parameterized by
    %   either mean parameters or information parameters
    %
    %   Mean parameters include
    %   - mu:       the mean vector 
    %   - C:        the covariance matrix object
    %
    %   Information parameters include
    %   - c0:      the constant term
    %   - h:       the potential vector
    %   - J:       the information matrix
    %  
    %   These two types of parameterization are related to each other
    %   as follows:
    %   - J = inv(sigma)
    %   - h = inv(sigma) * mu
    %   - c0 = mu' * inv(sigma) * mu 
    %   Or, equivalently,
    %   - sigma = inv(J)
    %   - mu = inv(J) * h
    %
    %   For Gaussian distributions with zero mean, one can simply set
    %   mu or h to a scalar 0, despite the actual dimension of 
    %   the space.
    %   
    %   Generally, the evaluation of Mahalanobis distance or probability
    %   density function, and the posterior computation relies on the
    %   information parameters, while sampling relies on mean parameters.
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on June 19, 2010
    %       - Modified by Dahua Lin, on Sep 15, 2010
    %       - Modified by Dahua Lin, on Nov 12, 2010
    %       - Modified by Dahua Lin, on Aug 14, 2011    
    %       - Modified by Dahua Lin, on Aug 25, 2011
    %
    
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')        
        dim;    % the dimension of the vector space (d)
        num;    % the number of models contained in the object (n)
                                   
        mu;     % the mean vector(s)
                % in general, it is a d x n matrix
                % when there is one model with zero mean vector, 
                % mu can be represented as a zero scalar
                
        C;      % the covariance matrix(ces) in a possibly compressed form
                % cform == 's': scalar or 1 x n row vector
                % cform == 'd': d x 1 vector or d x n matrix
                % cform == 'f': d x d matrix or d x d x n array
        
        c0;     % the constant term(s) (= mu * J * mu) 
        h;      % the potential vector(s) (= J * mu)
        J;      % the information matrix(ces) (= inv(cov))
                % Note: J can be represented in a compressed form 
                % depending on the value of cform
                
        ldcov;  % the value of log(det(cov))
        
        has_mp = false; % whether the mean parameters are available
        has_ip = false; % whether the information parameters are available
        shared_cov = false; % whether the covariance is shared
        zmean = false;  % test whether it is a single model with zero mean
    end
    
    
    %% constructor    
    
    methods(Static)
        
        function G = from_mp(mu, C, use_ip)
            % Create Gaussian distribution(s) from mean parameterization
            %
            %   G = gaussd.from_mp(mu, C);
            %       creates an object to represent Gaussian distributions
            %       using mean parameterization.
            %
            %       Input:
            %       - cf:   the form of covariance matrix
            %       - mu:   the mean vector(s) [d x n matrix]. If the
            %               mu is a zero vector, then it can be given
            %               as a zero scalar.
            %       - C:    the covariance matrix represented by a
            %               pdmat struct.            
            %
            %   G = gaussd.from_mp(mu, C, 'ip');
            %       additionally derives the information parameters.
            %                                    
            
            % verify input types
            
            if ~(isfloat(mu) && isreal(mu) && ndims(mu) == 2 && ~isempty(mu))
                error('gaussd:from_mp:invalidarg', ...
                    'mu should be a real matrix.');
            end
            
            if ~is_pdmat(C)
                error('gaussd:from_mp:invalidarg', ...
                    'C should be a pdmat struct.');
            end
            
            % determine & verify d and n
            
            d = C.d;
            
            if isequal(mu, 0)
                n = 1;                
            else
                n = size(mu, 2);
                if size(mu, 1) ~= d
                    error('gaussd:invalidarg', ...
                        'The dimensions of mu and C are inconsistent.');
                end
                if n == 1 && all(mu == 0)
                    mu = 0;
                end
            end
            
            if ~(C.n == n || C.n == 1)
                error('gaussd:invalidarg', ...
                    'C.n is not consistent with the actual number of models.');
            end
            
            sn = C.n;
            
            % determine whether to use ip (information parameter)
            
            if nargin >= 3
                if strcmpi(use_ip, 'ip')
                    uip = true;
                else
                    error('gaussd:from_mp:invalidarg', ...
                        'The 3rd argument can only be ''ip''.');
                end
            else
                uip = false;
            end
            
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;  
            
            G.mu = mu;                
            G.C = C;
            G.has_mp = true;
            
            G.shared_cov = (sn == 1);
            G.zmean = isequal(mu, 0);
                            
            % derive canonical parameters (if requested)
            
            if uip                
                Jm = pdmat_inv(C);
                ldc = pdmat_lndet(C);
                
                if G.zmean
                    hv = 0;
                    cv = 0;
                else
                    hv = pdmat_mvmul(Jm, mu);                    
                    cv = dot(hv, mu, 1);
                end
                
                G.c0 = cv;
                G.h = hv;
                G.J = Jm;
                G.ldcov = ldc;
                G.has_ip = true;                
            end
        end        
        
        
        function G = from_ip(h, J, c0, use_mp)
            % Create Gaussian distribution(s) from information parameters
            %
            %   G = gaussd.from_ip(h, J);
            %   G = gaussd.from_ip(h, J, c0);
            %       creates an object to represent Gaussian distributions
            %       using information parameterization.
            %
            %       Input:
            %       - h:    the potential vector(s) [d x n matrix]
            %       - J:    the information matrix [a pdmat struct]
            %       - c0:   the pre-computed constant term. 
            %               The function will compute it if not given.
            %
            %   G = gaussd.from_ip(h, J, [], 'mp');
            %   G = gaussd.from_ip(h, J, c0, 'mp');
            %       additionally derives the mean parameters.
            %
            
            % verify input types
            
            if ~(isfloat(h) && isreal(h) && ndims(h) == 2)
                error('gaussd:from_ip:invalidarg', ...
                    'h should be a real matrix.');
            end
            
            if ~(is_pdmat(J))
                error('gaussd:from_ip:invalidarg', ...
                    'J should be a pdmat struct');
            end
            
            if nargin < 3
                c0 = [];
            end
            
            if nargin >= 4
                if strcmpi(use_mp, 'mp')
                    ump = true;
                else
                    error('gaussd:from_ip:invalidarg', ...
                        'The 4th argument can only be ''mp''.');
                end
            else
                ump = false;
            end
            
            % determine d and n
            
            d = J.d;
            
            if isequal(h, 0)
                n = 1;                
            else
                n = size(h, 2);
                if size(h, 1) ~= d
                    error('gaussd:invalidarg', ...
                        'The dimensions of h and J are inconsistent.');
                end
                if n == 1 && all(h == 0)
                    h = 0;
                end
            end
            
            if ~(J.n == n || J.n == 1)
                error('gaussd:invalidarg', ...
                    'J.n is not consistent with the actual number of models.');
            end
            
            sn = J.n;
                        
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;
                    
            G.h = h;
            G.J = J;
            G.ldcov = - pdmat_lndet(J);
            G.has_ip = true;
            
            G.shared_cov = (sn == 1);
            G.zmean = isequal(h, 0);
            
            % derive mean parameters                     
            
            if ump                
                Cm = pdmat_inv(J);
                if G.zmean
                    mv = 0;                    
                else
                    mv = pdmat_mvmul(Cm, h);
                end 
            end
                        
            % determine c0 and ldc
                                                                                    
            if ~isempty(c0)
                if ~(isfloat(c0) && isreal(c0) && isequal(size(c0), [1 n]))
                    error('gaussd:from_ip:invalidarg', ...
                        'c0 should be a real vector of size 1 x n.');
                end
            else                
                if G.zmean
                    c0 = 0;
                else
                    mv = pdmat_lsolve(J, h);
                    c0 = dot(h, mv, 1);                    
                end
            end          
            G.c0 = c0;
            
            % construct Gaussian object         
            
            if ump
                G.mu = mv;
                G.C = Cm;
                G.has_mp = true;
            end           
        end
    end
    
    
    %% Statistics
    
    methods
        
        function V = mean(G)
            % Get the mean vectors of the distributions
            %
            %   V = mean(G);
            %       returns the mean vectors as columns of V.
            %       The size of V is [dim, num]
            %
            
            if ~G.has_mp
                error('gaussd:nomp', 'Mean parameterization is needed.');
            end
            
            if G.zmean
                V = zeros(G.dim, G.num);
            else
                V = G.mu;
            end
        end
        
        function V = var(G)
            % Get the marginal variances of all components
            %
            %   V = mean(G);
            %       returns the marginal variances as columns of V.
            %       The size of V is [dim, num]
            %
            
            if ~G.has_mp
                error('gaussd:nomp', 'Mean parameterization is needed.');
            end
            
            V = pdmat_diag(G.C);
            if G.shared_cov && G.num > 1
                V = repmat(V, [1, G.num]);
            end
        end
        
        
        function Cmat = cov(G, i)
            % Gets the covariance matrix of a distribution in G
            %
            %   Cmat = cov(G);
            %       returns the covariance matrix of G. This applies
            %       only to the case where G has only one covariance,
            %       or the covariance is shared.
            %
            %   Cmat = cov(G, i);
            %       returns the covariance matrix associated with the
            %       i-th distribution.
            %
            
            if ~G.has_mp
                error('gaussd:nomp', 'Mean parameterization is needed.');
            end            
            
            if nargin < 2
                if ~G.shared_cov
                    error('gaussd:invalidarg', ...
                        'Index needs to be given when there are multi-cov.');
                end
                
                Cmat = pdmat_fullform(G.C);
            else
                if G.shared_cov
                    Cmat = pdmat_fullform(G.C);
                else
                    Cmat = pdmat_fullform(G.C, i);
                end
            end                            
        end    
        
        
        function v = entropy(G, i)
            % Computes the entropy
            %
            %   v = entropy(G);
            %       computes the entropies of all distributions.
            %
            %   v = entropy(G, i);
            %       computes the entropies of the selected distributions.
            %               
            
            if G.has_mp
                is_inv = 0;
                a = G.C;
            else
                is_inv = 1;
                a = G.J;
            end
            
            if nargin < 2
                n = G.num;
                if n == 1
                    if is_inv
                        v = gentropy(a, 'inv');
                    else
                        v = gentropy(a);
                    end
                else
                    if is_inv
                        v = gentropy(a, 'inv');
                    else
                        v = gentropy(a);
                    end
                    if G.shared_cov
                        v = v(1, ones(1, n));
                    end
                end
            else
                if ~(isnumeric(i) && isvector(i))
                    error('gaussd:invalidarg', 'i should be a numeric vector.');
                end                
                if G.shared_cov
                    if is_inv
                        v = gentropy(a, 'inv');
                    else
                        v = gentropy(a);
                    end
                    n = numel(i);
                    if n > 1
                        v = v(1, ones(1, n));
                    end
                else
                    if is_inv
                        v = gentropy(pdmat_pick(a, i), 'inv');
                    else
                        v = gentropy(pdmat_pick(a, i));
                    end
                end
            end
        end
    
    end
    
    
    %% Probability evaluation
    
    methods
        
        function D = sqmahdist(G, X, si)
            % Compute the squared Mahalanobis distances from samples to centers
            %
            %   D = sqmahdist(G, X);
            %       computes the squared Mahalanobis distances between the 
            %       samples in X and the Gaussian centers (with respect
            %       to the corresponding covariance matrices).
            %
            %   D = sqmahdist(G, X, si);
            %       computes the squared Mahalanobis distances between the
            %       samples in X and the centers of the Gaussian 
            %       distributions selected by indices si.
            %
            %   Note that information parameters are required for the 
            %   computation.
            %
            
            if ~G.has_ip
                error('gaussd:sqmahdist:noip', ...
                    'Information parameters are needed.');
            end
            
            % take parameters
            
            c0v = G.c0;
            hv = G.h;
            Jm = G.J;
            zm = G.zmean;
            scov = G.shared_cov;
            
            if nargin >= 3 && ~isempty(si)
                c0v = c0v(1, si);
                if ~zm
                    hv = hv(:, si);
                end
                if ~scov
                    Jm = pdmat_pick(Jm, si);
                end
            end                        
                
            % compute 
            
            t2 = pdmat_quad(Jm, X, X);
            
            if zm
                D = t2;
            else
                t1 = hv' * X;
                
                n1 = size(t1, 1);
                n2 = size(t2, 1);
                if n1 == 1
                    D = t2 - 2 * t1 + c0v;
                else
                    if n2 == 1
                        D = bsxfun(@minus, t2, 2 * t1);
                    else
                        D = t2 - 2 * t1;
                    end
                    D = bsxfun(@plus, D, c0v.');
                end
            end
        end
        
        
        function L = logpdf(G, X, si)
            % Compute logarithm of PDF of given samples
            %
            %   L = logpdf(G, X)
            %       compute logarithm of probability density function
            %       at the samples given by columns of X.
            %
            %       Let m be the number of distributions contained in G,
            %       and n be the number of samples in X. Then L will be
            %       a matrix of size m x n, with L(i, j) being the pdf
            %       value at X(:,j) w.r.t the j-th Gaussian.
            %
            %   L = logpdf(G, X, si);
            %       compute the logarithm of pdf with respect to the 
            %       models selected by the index vector si.
            %
            %   Note that canonical parameters are required for the 
            %   computation.
            %
            
            if nargin < 3 || isempty(si)
                D = sqmahdist(G, X);
                a0 = G.ldcov + G.dim * log(2 * pi);
            else
                D = sqmahdist(G, X, si);
                if G.shared_cov
                    ldc = G.ldcov;
                else
                    ldc = G.ldcov(1, si);
                end
                a0 = ldc + G.dim * log(2 * pi);
            end
            
            if isscalar(a0)
                L = -0.5 * (D + a0);
            else
                L = -0.5 * bsxfun(@plus, D, a0.');
            end
        end
        
                                     
        function P = pdf(G, X, si)
            % Compute the probability density function
            %                        
            %   P = pdf(G, X)
            %       compute probability density function at the samples 
            %       given by columns of X.
            %
            %       Let m be the number of distributions contained in G,
            %       and n be the number of samples in X. Then L will be
            %       a matrix of size m x n, with P(i, j) being the pdf
            %       value at X(:,j) w.r.t the j-th Gaussian.
            %
            %   P = pdf(G, X, si);
            %       compute the probability density with respect to the 
            %       models selected by the index vector si.
            %
            %   Note that canonical parameters are required for the 
            %   computation.  
            %
            
            if nargin < 3
                P = exp(logpdf(G, X));
            else
                P = exp(logpdf(G, X, si));
            end            
        end                              
                
    end
    
    
    %% Inference and Sampling
    
    methods
        function [hp, Jp] = add_info(G, ho, Jo, i)
            % Adds information params to the current model
            %
            %   [hp, Jp] = G.add_info(ho, Jo);
            %       adds the information parameters that reflect the 
            %       feedback of the observations, to obtain the posterior
            %       information parameters.
            %
            %   [hp, Jp] = G.add_info(ho, Jo, i);
            %       the information are to be added to the i-th model
            %       in this object.
            %
            
            % verify input arguments
            
            if ~G.has_ip
                error('gaussd:add_info:noip', ...
                    'Information parameters are needed.');
            end
            
            [d, n] = size(ho);
            if d ~= G.dim
                error('gaussd:add_info:invalidarg', ...
                    'The input dimension does not match the model dimension.');
            end                        
            if ~(Jo.n == 1 || Jo.n == n)
                error('gaussd:add_info:invalidarg', ...
                    'The #components of the inputs ho and Jo are inconsistent.');
            end 
                                                
            if nargin < 4 || isempty(i)
                n0 = G.num;
                if ~(n0 == 1 || n0 == n)
                    error('gaussd:add_info:invalidarg', ...
                        'The #components of the inputs do not match the models.');
                end
                                
                h0 = G.h;
                J0 = G.J;
            else
                if ~(isnumeric(i) && isscalar(i) && i == fix(i) && i >= 1)
                    error('gaussd:add_info:invalidarg', ...
                        'The index i should be an integer scalar.');
                end
                if n ~= 1
                    error('gaussd:add_info:invalidarg', ...
                        'The #components of input must be one when i is specified.');
                end
                
                n0 = 1;
                h0 = G.h(:, i);
                if G.shared_cov
                    J0 = G.J;
                else
                    J0 = pdmat_pick(G.J, i);
                end
            end
             
            % main
            
            if n == n0
                hp = h0 + ho;            
            else
                hp = bsxfun(@plus, h0, ho);
            end
            
            Jp = pdmat_plus(J0, cf0, Jo, cfo); 
        end
        
        
        function Mp = pos_mean(G, ho, Jo, i)
            % Get the posterior mean(s) given observed information
            %
            %   Mp = G.pos_mean(ho, Jo);
            %   Mp = G.pos_mean(ho, Jo, i);
            %
            
            if nargin < 5
                [hp, Jp] = add_info(G, ho, Jo);
            else
                [hp, Jp] = add_info(G, ho, Jo, i);
            end
            
            Mp = pdmat_lsolve(Jp, hp);            
        end
        
        
        function X = pos_sample(G, ho, Jo, n, i)
            % sample from posterior Gaussian distribution
            %
            %   X = pos_sample(G, ho, Jo, n);
            %   X = pos_sample(G, ho, Jo, n, i);
            %       draws n samples from the model given observed
            %       information.
            %
            
            if size(ho, 2) ~= 1
                error('gaussd:pos_sample:invalidarg', ...
                    'The #components of the input must be 1.');
            end
            
            if nargin < 5
                if G.num > 1
                    error('gaussd:pos_sample:invalidarg', ...
                        'i is needed when G contains multiple models.');
                end
                [hp, Jp] = add_info(G, ho, Jo);
            else
                if ~(isnumeric(i) && isscalar(i) && i == fix(i) && i >= 1)
                    error('gaussd:pos_sample:invalidarg', ...
                        'i must be a positive integer scalar.');
                end
                [hp, Jp] = add_info(G, ho, Jo, i);
            end
            
            C_p = pdmat_inv(Jp);
            mu_p = pdmat_mvmul(C_p, hp);
            
            X = gsample(mu_p, C_p, n);            
        end        
        
        
        function X = sample(G, n, i)
            % Samples from the Gaussian distribution
            %            
            %   X = G.sample();
            %   X = G.sample(n);
            %       draws n samples from the Gaussian distribution. 
            %       (It must be G.num == 1 for this syntax).
            %
            %       When n is omitted, it is assumed to be 1.
            %
            %   X = G.sample(n, i);
            %       draws n samples from the i-th distribution.
            %
            %       When i is a vector, then n can be a vector of the 
            %       same size. In this case, it draws n(j) samples from
            %       the i(j)-th distribution.                     
            %
            
            if ~G.has_mp
                error('gaussd:sample:nomp', ...
                    'Mean parameters are needed.');
            end            
            if nargin < 2; n = 1; end
            
            if nargin < 3 || isempty(i)
                if G.num > 1
                    error('gaussd:sample:invalidarg', ...
                        'i is needed when G contains multiple models.');
                end                
                X = gsample(G.mu, G.C, n);
            else
                if ~(isvector(n) && isnumeric(n))
                    error('gaussd:pos_sample:invalidarg', ...
                        'n must be a numeric vector');
                end
                if ~(isvector(i) && isnumeric(i))
                    error('gaussd:pos_sample:invalidarg', ...
                        'i must be a numeric vector');
                end
                if numel(n) ~= numel(i)
                    error('gaussd:pos_sample:invalidarg', ...
                        'The sizes of n and i are inconsistent.');
                end               
                
                mu_ = G.mu;
                C_ = G.C;
                scov = G.shared_cov;
                
                if isscalar(i)
                    if scov
                        X = gsample(mu_, C_, n);
                    else
                        X = gsample(mu_(:,i), pdmat_pick(C_, i), n);
                    end
                else
                    if scov                    
                        X1 = gsample(mu_(:,i(1)), C_, n(1));
                    else
                        X1 = gsample(mu_(:,i(1)), pdmat_pick(C_, i(1)), n(1));
                    end
                    N = sum(n);
                    X = zeros(G.dim, N, class(X1));
                    X(:, 1:n(1)) = X1;
                    
                    ek = n(1);
                    for k = 2 : numel(n)
                        if scov
                            Xk = gsample(mu_(:,i(k)), C_, n(k));
                        else
                            Xk = gsample(mu_(:,i(k)), pdmat_pick(C_, i(k)), n(k));
                        end
                        
                        sk = ek + 1;
                        ek = ek + n(k);
                        X(:, sk:ek) = Xk;
                    end
                end
            end
        end
                
    end
    
    
    methods
        
        %% Visualization
        
        function plot_ellipse(G, r, varargin)
            % Plot ellipse representing the Gaussian models
            %            
            %   plot_ellipse(G, r, ...);
            %
            %       This function only works when dim == 2.
            %       
            
            if ~(G.has_mp && G.dim == 2)
                error('gaussd:plot_ellipse:invalidarg', ...
                    'The model should have G.has_mp and G.dim == 2.');
            end            
            
            mu_ = double(G.mu);
            C_ = double(G.C); 
            n = G.num;
            scov = G.shared_cov;
            
            for i = 1 : n
                
                if scov || n == 1
                    cc = C_;
                else
                    cc = pdmat_sub(C_, i);
                end                
                u = mu_(:, i);                                
                
                ns = 500;
                t = linspace(0, 2*pi, ns);
                x0 = [cos(t); sin(t)];
                x = bsxfun(@plus, pdmat_choltrans(cc, x0) * r, u);
                
                hold on;
                plot(x(1,:), x(2,:), varargin{:});                
            end
            
        end
        
        
    end
        
end



