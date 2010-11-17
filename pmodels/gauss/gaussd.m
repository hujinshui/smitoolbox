classdef gaussd
    % The class to represent Gaussian distribution(s)
    %
    %   Each object of this class can contain one or multiple Gaussian
    %   distributions. A Gaussian distribution can be parameterized by
    %   either mean parameters or information parameters
    %
    %   Mean parameters include
    %   - mu:       the mean vector 
    %   - sigma:    the covariance matrix object
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
    %           
    %
    
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')        
        dim;    % the dimension of the vector space
        num;    % the number of models contained in the object
                        
        mu;     % the mean vector(s)
        sigma;  % the covariance matrix object(s)
        
        c0;     % the constant term(s) in information parameters
        h;      % the potential vector
        J;      % the information matrix       
        ldcov;  % the value of log(det(cov))
        
        has_mp = false; % whether the mean parameters are available
        has_ip = false; % whether the information parameters are available
    end
    
    
    %% constructor    
    
    methods(Static)
        
        function G = from_mp(mu, sigma, uip)
            % Create Gaussian distribution(s) from mean parameterization
            %
            %   G = gaussd.from_mp(mu, sigma);
            %       creates an object to represent Gaussian distributions
            %       using mean parameterization.
            %
            %       Input:
            %       - mu:   the mean vector(s) [d x n matrix]
            %               if the mean vector is zero, if can be simply
            %               specified as a scalar 0 (in this case, sigma
            %               can contain only one matrix).
            %       - sigma: the covariance matrix. It should be an
            %                object of a symmetric matrix class.
            %                (udmat, dmat, gsymat, etc).
            %
            %   G = gaussd.from_mp(mu, sigma, 'ip');
            %       additionally derives the information parameters.
            %                                    
            
            % verify input types
            
            if ~(isfloat(mu) && ndims(mu) == 2)
                error('gaussd:from_mp:invalidarg', ...
                    'mu should be a numeric matrix.');
            end
            
            if ~(isobject(sigma))
                error('gaussd:from_mp:invalidarg', ...
                    'sigma should be an object of symmetric matrix.');
            end
            
            % determine d and n
            
            if isequal(mu, 0)
                d = sigma.d;
                if sigma.n ~= 1
                    error('gaussd:from_mp:invalidarg', ...
                        'sigma must be a single-matrix object when mu is 0.');
                end
                n = 1;
            else
                [d, n] = size(mu);
                if sigma.d ~= d
                    error('gaussd:from_mp:invalidarg', ...
                        'The dimension of sigma does not match that of mu.');
                end
                
                sn = sigma.n;
                if sn > 1 && sn ~= n
                    error('gaussd:from_mp:invalidarg', ...
                        'The sigma.n does not match size(mu, 2).');
                end
            end
            
            % determine whether to use ip
            
            if nargin >= 3
                if strcmp(uip, 'ip')
                    uip = true;
                else
                    error('gaussd:from_mp:invalidarg', ...
                        'The 3rd argument should be ''ip''.');
                end
            else
                uip = false;
            end
            
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;
            G.mu = mu;
            G.sigma = sigma;
            G.has_mp = true;
                            
            % derive canonical parameters (if requested)
            
            if uip
                Jm = inv(sigma);
                ldc = lndet(sigma);
                
                if isequal(mu, 0)
                    hv = 0;
                    cv = 0;
                else
                    if sn == 1
                        hv = Jm * mu; %#ok<MINV>
                    else
                        hv = cmv(Jm, mu);
                    end
                    cv = dot(hv, mu, 1);
                end
                
                G.c0 = cv;
                G.h = hv;
                G.J = Jm;
                G.ldcov = ldc;
                G.has_ip = true;                
            end
        end        
        
        
        function G = from_ip(h, J, c0, ump)
            % Create Gaussian distribution(s) from information parameters
            %
            %   G = gaussd.from_ip(h, J);
            %   G = gaussd.from_ip(h, J, c0);
            %       creates an object to represent Gaussian distributions
            %       using information parameterization.
            %
            %       Input:
            %       - h:    the potential vectors [d x n matrix]
            %               For the models with zero means, h can be
            %               simply input as a scalar 0.            
            %       - J:    the information matrix object
            %       - c0:   the constant term. The function will compute
            %               c0 if it is not specified.
            %
            %   G = gaussd.from_ip(h, J, [], 'mp');
            %   G = gaussd.from_ip(h, J, c0, 'mp');
            %       additionally derives the mean parameters.
            %
            
            % verify input types
            
            if ~(isfloat(h) && ndims(h) == 2)
                error('gaussd:from_cp:invalidarg', ...
                    'h should be a numeric matrix.');
            end
            
            if ~isobject(J)
                error('gaussd:from_cp:invalidarg', ...
                    'J should be an object of symmetric matrix.');
            end
            
            if nargin < 3
                c0 = [];
            end
            
            if nargin >= 4
                if strcmp(ump, 'mp')
                    ump = true;
                else
                    error('gaussd:from_cp:invalidarg', ...
                        'The 4th argument should be ''mp''.');
                end
            else
                ump = false;
            end
            
            % determine d and n
            
            if isequal(h, 0)
                d = J.d;
                n = J.n;             
            else
                [d, n] = size(h);
                if J.d ~= d
                    error('gaussd:from_cp:invalidarg', ...
                        'The dimension of J does not match that of h.');
                end
                
                sn = J.n;
                if sn > 1 && sn ~= n
                    error('gaussd:from_cp:invalidarg', ...
                        'The J.n does not match size(h, 2).');
                end
            end
            
            
            % determine mu and sigma                       
            
            if ump                
                sigma_ = inv(J);
                if isequal(h, 0)
                    mu_ = 0;
                else
                    if sn == 1
                        mu_ = sigma_ * h; %#ok<MINV>
                    else
                        mu_ = cmv(sigma_, h);
                    end
                end                
            else
                if isequal(h, 0)
                    mu_ = 0;
                else
                    if sn == 1
                        mu_ = J \ h;
                    else
                        mu_ = cdv(J, h);
                    end
                end
            end
                        
            % determine c0 and ldc
                                                                                    
            if ~isempty(c0)
                if ~(isfloat(c0) && isequal(size(c0), [1 n]))
                    error('gaussd:from_cp:invalidarg', ...
                        'c0 should be a numeric vector of size 1 x n.');
                end
            else                
                if isequal(h, 0)
                    c0 = 0;
                else
                    c0 = dot(h, mu_, 1);
                end
            end                      
            
            ldc = -lndet(J);
            
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;
            G.c0 = c0;
            G.h = h;
            G.J = J;
            G.ldcov = ldc;
            
            G.has_ip = true;
            
            if ump
                G.mu = mu_;
                G.sigma = sigma_;
                G.has_mp = true;
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
                error('gaussd:sqmahdist:nocp', ...
                    'Information parameters are required.');
            end
            
            % take parameters
            
            c0_ = G.c0;
            h_ = G.h;
            J_ = G.J;
            
            if nargin >= 3
                c0_ = c0_(:, si);
                if ~isequal(h_, 0)
                    h_ = h_(:, si);
                end
                if J_.n ~= 1
                    J_ = J_.take(si);
                end
            end                        
                
            % compute 
            
            t2 = quad(J_, X, X);
            
            if isequal(h_, 0)
                D = t2;
            else
                if isequal(h_, 0)
                    t1 = 0;
                else
                    t1 = h_' * X;
                end
                
                n1 = size(h_, 2);
                n2 = J_.n;
                if n1 == 1
                    D = t2 - 2 * t1 + c0_;
                else
                    if n2 == 1
                        D = bsxfun(@minus, t2, 2 * t1);
                    else
                        D = t2 - 2 * t1;
                    end
                    D = bsxfun(@plus, D, c0_.');
                end
            end
        end
        
        
        function D = sqmahdist_map(G, X, M)
            % Compute the squared Mahalanobis distances to remapped centers
            %
            %   D = G.sqmahdist_remap(X, M);
            %       
            %       Suppose X comprises n samples, then M should be
            %       a vector of 1 x n, and the output D is also a
            %       vector of size 1 x n.
            %
            %       D(i) equals the squared Mahalanobis distance between
            %       X(:,i) the the M(i)-th center of the object.
            %
            
            if ~G.has_ip
                error('gaussd:sqmahdist_remap:nocp', ...
                    'Information parameters are required.');
            end
            
            nx = size(X, 2);
            if nx ~= size(M, 2)
                error('gaussd:sqmahdist_remap:invalidarg', ...
                    'The sizes of X and M are inconsistent.');
            end            
            
            % take parameters
            
            c0_ = G.c0;
            h_ = G.h;
            J_ = G.J;
            
            % compute
            
            if isequal(h_, 0)
                t1 = 0;
            else
                t1 = dot(h_(:, M), X, 1);
            end
            
            if J_.n == 1
                t2 = quad(J_, X, X);                
            else                
                K = J_.n;
                t2 = zeros(1, nx);
                
                if nx <= 2 * K
                    for i = 1 : nx
                        cx = X(:, i);
                        t2(i) = quad(J_.take(M(i)), cx, cx);
                    end
                else
                    gs = intgroup(K, M);                
                    for k = 1 : K
                        cg = gs{k};
                        cX = X(:,cg);
                        t2(cg) = quad(J_.take(k), cX, cX);
                    end
                end
            end
            
            if ~isscalar(c0_)
                c0_ = c0_(M);
            end
            
            D = t2 - 2 * t1 + c0_;
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
            
            if nargin < 3             
                D = sqmahdist(G, X);
                a0 = G.ldcov + G.dim * log(2 * pi);
            else
                D = sqmahdist(G, X, si);
                if G.J.n == 1
                    a0 = G.ldcov + G.dim * log(2 * pi);
                else
                    a0 = G.ldcov(:, si) + G.dim * log(2 * pi);
                end
            end
            
            if isscalar(a0)
                L = -0.5 * (D + a0);
            else
                L = -0.5 * bsxfun(@plus, D, a0.');
            end
        end
        
        
        function L = logpdf_map(G, X, M)
            % Compute log-pdf of samples with respect to mapped Gaussians
            %
            %   L = logpdf_map(G, X, M);
            %       
            %       Suppose X comprises n samples, then M should be
            %       a vector of 1 x n, and the output D is also a
            %       vector of size 1 x n.
            %
            %       L(i) equals the log-pdf of X(:,i) with respect to
            %       the M(i)-the Gaussian contained in G.
            %
            
            D = sqmahdist_map(G, X, M);
            if G.J.n == 1
                a0 = G.ldcov + G.dim * log(2 * pi);
            else
                a0 = G.ldcov(M) + G.dim * log(2 * pi);
            end
            
            L = -0.5 * (D + a0);
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
        
        
        function P = pdf_map(G, X, M)
            % Compute pdf of samples with respect to mapped Gaussians
            %
            %   P = logpdf_map(G, X, M);
            %       
            %       Suppose X comprises n samples, then M should be
            %       a vector of 1 x n, and the output D is also a
            %       vector of size 1 x n.
            %
            %       P(i) equals the pdf of X(:,i) with respect to
            %       the M(i)-the Gaussian contained in G.
            %
            
            P = logpdf_map(G, X, M);
        end
    end
    
    
    %% Inference and Sampling
    
    methods
        
        
        function [h1, J1] = inject(G, ha, Ja, i)
            % Compute the posterior with injected observations
            %
            %   [h1, J1] = inject(G, ha, Ja);
            %   [h1, J1] = inject(G, ha, Ja, i);
            %       computes the information parameters of the posterior 
            %       Gaussian distribution with the observations 
            %       injected with ha and Ja, which are quantities to be 
            %       added to h and J respectively.            
            %
            %       The output h1 and J1 are respectively the potential
            %       vector and information matrix of the posterior
            %       Gaussian distribution.
            %
            
            if ~G.has_ip
                error('gaussd:inject:nocp', 'Information parameters are required.');
            end
            
            h0 = G.h;
            J0 = G.J;
            if nargin >= 4 && ~isempty(i)
                if ~isequal(h0, 0)
                    h0 = h0(:, i);
                end
                
                if J0.n ~= 1
                    J0 = J0.take(i);
                end
            end
            
            if isequal(h0, 0)
                h1 = ha;
            else
                if size(h0, 2) == size(ha, 2)
                    h1 = h0 + ha;
                else
                    h1 = bsxfun(@plus, h0, ha);
                end
            end
            
            J1 = J0 + Ja; 
        end
        
        
        function Gpos = posterior(G, ha, Ja, i, ump)
            % Compute the posterior Gaussian distribution
            %
            %   Gpos = posterior(G, ha, Ja);
            %   Gpos = posterior(G, ha, Ja, i);
            %       Gets the posterior Gaussian given the observations
            %       summarized by ha and Ja, which are to be injected
            %       to the information parameters of the prior.
            %       
            %       By default, the returned object is with only
            %       information parameters.                        
            %   
            %   Gpos = posterior(G, ha, Ja, [], 'mp');
            %   Gpos = posterior(G, ha, Ja, i, 'mp');                  
            %       Returns the posterior distribution object with
            %       both information parameters together with
            %       mean parameters.
            %
            
            if nargin < 4 || isempty(i)
                [h1, J1] = inject(G, ha, Ja);
            else
                [h1, J1] = inject(G, ha, Ja, i);
            end
            
            if nargin < 5
                Gpos = gaussd.from_ip(h1, J1);
            else
                Gpos = gaussd.from_ip(h1, J1, [], ump);
            end                        
        end
        
        
        function [M, C] = pos_mean(G, ha, Ja, i)
            % Compute the mean parameters of posterior distribution
            %
            %   M = G.posterior_mean(G, ha, Ja);
            %   M = G.posterior_mean(G, ha, Ja, i);
            %
            %       solves the mean of the posterior Gaussian distribution.
            %       
            %   [M, C] = G.posterior_mean(G, ha, Ja);
            %   [M, C] = G.posterior_mean(G, ha, Ja, i);
            %
            %       additionally returns the covariance object.
            %
            
            if nargin < 4 || isempty(i)
                [h1, J1] = inject(G, ha, Ja);
            else
                [h1, J1] = inject(G, ha, Ja, i);
            end
            
            if nargout < 2
                M = cdv(J1, h1);
            else
                C = inv(J1);
                M = cmv(C, h1);
            end            
        end

        
        function X = sample(G, n, i, rstream)
            % Samples from the Gaussian distribution
            %
            %   X = G.sample();                                
            %       draws a sample from the Gaussian distribution G.
            %       
            %       If G contains m distributions, then it draws
            %       m samples in total, one from each distribution.
            %       In particular, X will be a d x m matrix, and
            %       X(:,k) is from the k-th distribution.
            %
            %   X = G.sample(n);
            %       draws n samples from each distribution. 
            %   
            %       If there are m distributions in G, then the columns 
            %       in X(:,1:m) are from the 1st distribution, and
            %       X(:, m+1:2*m) are from the 2nd, and so on.            
            %
            %   X = G.sample(n, i);
            %       draws n samples from the i-th distribution.
            %
            %       Here, i can be a vector. For example, if i is [2, 3], 
            %       then the function draws n samples from the 2nd 
            %       distribution, and then draws n samples from the
            %       3rd one.
            %
            %       When i is a vector, then n can be a vector of the 
            %       same size. In this case, it draws n(j) samples from
            %       the i(j)-th distribution.
            %
            %   X = G.sample(n, [], rstream);
            %   X = G.sample(n, i, rstream);
            %       one can further specifies the random number stream
            %       to be used in the sampling.            
            %
            
            if ~G.has_mp
                error('gaussd:sample:invalidarg', ...
                    'Mean parameters are required.');
            end
            
            if nargin < 2; n = 1; end
            if nargin < 3; i = []; end
            if nargin < 4; rstream = []; end
            
            X = gsample(G.mu, G.sigma, n, i, rstream);         
        end
        
        
        function X = pos_sample(G, ha, Ja, n, i, rstream)
            % Samples from posterior distribution                  
            %
            %   X = G.pos_sample(ha, Ja);
            %   X = G.pos_sample(ha, Ja, n, i);             
            %   X = G.pos_sample(ha, Ja, n, [], rstream);
            %   X = G.pos_sample(ha, Ja, n, i, rstream);
            %
            %       draws a sample from the posterior distribution 
            %       conditioned on the observations injected through
            %       ha and Ja (addends to the potential vector and
            %       information matrix).                        
            %
            %       please refer to the help of the sample method for 
            %       details of what different syntax means.
            %

            if nargin < 4; n = 1; end
            if nargin < 5; i = []; end
            if nargin < 6; rstream = []; end
                        
            [M, C] = pos_mean(G, ha, Ja, i);
            X = gsample(M, C, n, [], rstream);
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
            
            mu_ = G.mu;
            sigma_ = G.sigma;
            
            for i = 1 : G.num
                
                if sigma_.n == 1
                    C = sigma_;
                else
                    C = sigma_.take(i);
                end                
                u = mu_(:, i);                                
                
                ns = 500;
                t = linspace(0, 2*pi, ns);
                x = bsxfun(@plus, C.choltrans([cos(t); sin(t)]) * r, u);
                
                hold on;
                plot(x(1,:), x(2,:), varargin{:});                
            end
            
        end
        
        
    end
        
end



