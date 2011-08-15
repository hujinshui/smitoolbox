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
    %
    
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')        
        dim;    % the dimension of the vector space (d)
        num;    % the number of models contained in the object (n)
        cform;  % the form of covariance matrix 
                % 's':  isotropic covariance represented by scalar(s)
                % 'd':  diagonal covariance represented by a vector of
                %       diagonal entries
                % 'f':  full covariance matrix(ces)
                                   
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
        
        function G = from_mp(cf, mu, C, use_ip)
            % Create Gaussian distribution(s) from mean parameterization
            %
            %   G = gaussd.from_mp(cf, mu, C);
            %       creates an object to represent Gaussian distributions
            %       using mean parameterization.
            %
            %       Input:
            %       - cf:   the form of covariance matrix
            %       - mu:   the mean vector(s) [d x n matrix]
            %       - C:    the covariance matrix represented in
            %               the form specified by cf
            %
            %       The form of C:
            %       cf == 's':  the covariance is like cov = c * I.
            %                   As input, C is a scalar or 1 x n row
            %                   vector, when n models are to be packed in
            %                   this object.
            %       cf == 'd':  the covariance is a diagonal matrix.
            %                   As input, C is a dx1 column vector to 
            %                   represent the diagonal entries, or 
            %                   a d x n matrix when there are n models.
            %       cf == 'f':  the covariance is in full matrix form.
            %                   As input, C is a d x d covariance matrix,
            %                   or a d x d x n array when with n models.
            %
            %   G = gaussd.from_mp(cf, mu, C, 'ip');
            %       additionally derives the information parameters.
            %                                    
            
            % verify input types
            
            if ~(ischar(cf) && isscalar(cf))
                error('gaussd:from_mp:invalidarg', ...
                    'cf should be a char scalar.');
            end
            
            if ~(isfloat(mu) && isreal(mu) && ndims(mu) == 2 && ~isempty(mu))
                error('gaussd:from_mp:invalidarg', ...
                    'mu should be a real matrix.');
            end
            
            if ~(isfloat(C) && isreal(C) && ~isempty(mu))
                error('gaussd:from_mp:invalidarg', ...
                    'C should be a real array.');
            end
            
            % determine d and n
            
            [d, n] = size(mu);
            
            % verify the size of C
            
            switch cf
                case 's'
                    csiz1 = [1, 1];
                    csiz = [1, n];
                case 'd'
                    csiz1 = [d, 1];
                    csiz = [d, n];
                case 'f'
                    csiz1 = [d, d];
                    if n == 1
                        csiz = [d, d];
                    else
                        csiz = [d, d, n];
                    end
                otherwise
                    error('gaussd:from_mp:invalidarg', ...
                        'The argument cf is invalid.');
            end
                    
            if isequal(size(C), csiz1)
                sn = 1;
            elseif isequal(size(C), csiz)
                sn = n;
            else
                error('gaussd:from_mp:invalidarg', ...
                    'The size of C is incorrect.');
            end                
            
            % determine whether to use ip (information parameter)
            
            if nargin >= 4
                if strcmpi(use_ip, 'ip')
                    uip = true;
                else
                    error('gaussd:from_mp:invalidarg', ...
                        'The 4th argument can only be ''ip''.');
                end
            else
                uip = false;
            end
            
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;
            G.cform = cf;   
            
            G.mu = mu;
            G.C = C;
            G.has_mp = true;
            
            G.shared_cov = (sn == 1);
            G.zmean = (n == 1 && all(mu == 0));
                            
            % derive canonical parameters (if requested)
            
            if uip                
                Jm = gmat_inv(cf, C);
                ldc = gmat_lndet(cf, d, C);
                
                if G.zmean
                    hv = 0;
                    cv = 0;
                else
                    hv = gmat_mvmul(cf, Jm, mu);                    
                    cv = dot(hv, mu, 1);
                end
                
                G.c0 = cv;
                G.h = hv;
                G.J = Jm;
                G.ldcov = ldc;
                G.has_ip = true;                
            end
        end        
        
        
        function G = from_ip(cf, h, J, c0, use_mp)
            % Create Gaussian distribution(s) from information parameters
            %
            %   G = gaussd.from_ip(cf, h, J);
            %   G = gaussd.from_ip(cf, h, J, c0);
            %       creates an object to represent Gaussian distributions
            %       using information parameterization.
            %
            %       Input:
            %       - cf:   the form of covariance matrix
            %       - h:    the potential vector(s) [d x n matrix]
            %       - J:    the information matrix in specified form (cf)
            %       - c0:   the pre-computed constant term. 
            %               The function will compute it if not given.
            %
            %   G = gaussd.from_ip(cf, h, J, [], 'mp');
            %   G = gaussd.from_ip(cf, h, J, c0, 'mp');
            %       additionally derives the mean parameters.
            %
            
            % verify input types
            
            if ~(ischar(cf) && isscalar(cf))
                error('gaussd:from_ip:invalidarg', ...
                    'cf should be a char scalar.');
            end
            
            if ~(isfloat(h) && isreal(h) && ndims(h) == 2)
                error('gaussd:from_ip:invalidarg', ...
                    'h should be a real matrix.');
            end
            
            if ~(isfloat(J) && isreal(J))
                error('gaussd:from_ip:invalidarg', ...
                    'J should be an array of real values.');
            end
            
            if nargin < 4
                c0 = [];
            end
            
            if nargin >= 5
                if strcmpi(use_mp, 'mp')
                    ump = true;
                else
                    error('gaussd:from_ip:invalidarg', ...
                        'The 5th argument can only be ''mp''.');
                end
            else
                ump = false;
            end
            
            % determine d and n
            
            [d, n] = size(h);
            
            % verify the size of J
            
            switch cf
                case 's'
                    jsiz1 = [1, 1];
                    jsiz = [1, n];
                case 'd'
                    jsiz1 = [d, 1];
                    jsiz = [d, n];
                case 'f'
                    jsiz1 = [d, d];
                    if n == 1
                        jsiz = [d, d];
                    else
                        jsiz = [d, d, n];
                    end
                otherwise
                    error('gaussd:from_ip:invalidarg', ...
                        'The argument cf is invalid.');
            end
                    
            if isequal(size(J), jsiz1)
                sn = 1;
            elseif isequal(size(J), jsiz)
                sn = n;
            else
                error('gaussd:from_ip:invalidarg', ...
                    'The size of J is incorrect.');
            end 
                        
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;
            G.cform = cf;
                    
            G.h = h;
            G.J = J;
            G.ldcov = - gmat_lndet(cf, d, J);
            G.has_ip = true;
            
            G.shared_cov = (sn == 1);
            G.zmean = (n == 1 && all(h == 0));
            
            % derive mean parameters                     
            
            if ump                
                Cm = gmat_inv(cf, J);
                if isequal(h, 0)
                    mv = 0;
                else
                    mv = gmat_mvmul(cf, Cm, h);
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
            
            c0_ = G.c0;
            h_ = G.h;
            J_ = G.J;
            zm = G.zmean;
            scov = G.shared_cov;
            
            if nargin >= 3
                c0_ = c0_(1, si);
                if ~zm
                    h_ = h_(:, si);
                end
                if ~scov
                    J_ = gmat_sub(cf, J_, si);
                end
            end                        
                
            % compute 
            
            t2 = gmat_quad(cf, J_, X, X);
            
            if zm
                D = t2;
            else
                t1 = h_' * X;
                
                n1 = size(t1, 1);
                n2 = size(t2, 1);
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
                error('gaussd:inject:noip', ...
                    'Information parameters are needed.');
            end
            
            h0 = G.h;
            J0 = G.J;
            zm = G.zmean;
            scov = G.shared_cov;
            
            if nargin >= 4 && ~isempty(i)
                if ~zm
                    h0 = h0(:, i);
                end
                
                if ~scov
                    J0 = J0.take(i);
                end
            end
            
            if zm
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
                Gpos = gaussd.from_ip(G.cf, h1, J1);
            else
                Gpos = gaussd.from_ip(G.cf, h1, J1, [], ump);
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
                M = gmat_lsolve(G.cf, J1, h1);
            else
                C = inv(J1);
                M = gmat_mvmul(G.cf, C, h1);
            end            
        end

        
        function X = sample(G, n, i)
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
            
            if ~G.has_mp
                error('gaussd:sample:nomp', ...
                    'Mean parameters are needed.');
            end
            
            if nargin < 2; n = 1; end
            if nargin < 3; i = []; end
            
            X = gsample(G.mu, G.C, n, i);         
        end
        
        
        function X = pos_sample(G, ha, Ja, n, i)
            % Samples from posterior distribution                  
            %
            %   X = G.pos_sample(ha, Ja);
            %   X = G.pos_sample(ha, Ja, n, i);             
            %
            %       draws a sample from the posterior distribution 
            %       conditioned on the observations injected through
            %       ha and Ja (addends to the potential vector and
            %       information matrix).                        
            %
            %

            if nargin < 4; n = 1; end
            if nargin < 5; i = []; end
                        
            [Mu, Cm] = pos_mean(G, ha, Ja, i);
            X = gsample(Mu, Cm, n, []);
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
            cf_ = G.cf;
            
            for i = 1 : n
                
                if scov || n == 1
                    cc = C_;
                else
                    cc = gmat_sub(cf_, C_, i);
                end                
                u = mu_(:, i);                                
                
                ns = 500;
                t = linspace(0, 2*pi, ns);
                x0 = [cos(t); sin(t)];
                x = bsxfun(@plus, gmat_choltrans(cf_, cc, x0) * r, u);
                
                hold on;
                plot(x(1,:), x(2,:), varargin{:});                
            end
            
        end
        
        
    end
        
end



