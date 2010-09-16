classdef gaussd
    % The class to represent Gaussian distribution
    %
    %   Each object of this class can contain one or multiple Gaussian
    %   distributions. A Gaussian distribution can be parameterized by
    %   either mean parameters or canonical parameters
    %
    %   Mean parameters are
    %   - mean:     the mean vector 
    %   - cov:      the covariance matrix object
    %
    %   Canonical parameterization is given by 
    %   - coef0:   the constant term
    %   - coef1:   the coefficient vector of first order term
    %   - coef2:   the coefficient matrix of second order term    
    %  
    %   These two types of parameterization are related to each other
    %   as follows:
    %   - coef2 = inv(cov)
    %   - coef1 = inv(cov) * mu
    %   - coef0 = mu' * inv(sigma) * mu 
    %   Or, equivalently,
    %   - cov = inv(coef2)
    %   - mean = cov * coef1
    %
    %   For Gaussian distributions with zero mean, one can simply set
    %   coef1 or mean to a scalar 0 despite the actual dimension of 
    %   the space.
    %   
    %   Generally, the evaluation of Mahalanobis distance or probability
    %   density function, and the posterior computation relies on the
    %   canonical parameters, while sampling relies on mean parameters.
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on June 19, 2010
    %       - Modified by Dahua Lin, on Sep 15, 2010
    %           
    %
    
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')        
        dim;    % the dimension of the vector space
        num;    % the number of models contained in the object
                        
        mean;   % the mean vector(s)
        cov;    % the covariance matrix object
        
        coef0;  % the constant term(s) in canonical parameters: mean' * inv(cov) * mean
        coef1;  % the coefficient vector(s) of the linear term: inv(cov) * mean
        coef2;  % the coefficient matrix object of the quadratic term: inv(cov)        
        ldcov;  % the value of log(det(cov))
        
        has_mp = false; % whether the mean parameters are available
        has_cp = false; % whether the canonical parameters are available
    end
    
    
    %% constructor    
    
    methods(Static)
        
        function G = from_mp(mu, sigma, ucp)
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
            %   G = gaussd.from_mp(mu, sigma, 'cp');
            %       additionally derives the canonical parameters.
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
            
            % determine whether to use cp
            
            if nargin >= 3
                if strcmp(ucp, 'cp')
                    ucp = true;
                else
                    error('gaussd:from_mp:invalidarg', ...
                        'The 3rd argument should be ''cp''.');
                end
            else
                ucp = false;
            end
            
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;
            G.mean = mu;
            G.cov = sigma;
            G.has_mp = true;
                            
            % derive canonical parameters (if requested)
            
            if ucp
                c2 = inv(sigma);
                ldc = lndet(sigma);
                
                if isequal(mu, 0)
                    c1 = 0;
                    c0 = 0;
                else
                    if sn == 1
                        c1 = c2 * mu; %#ok<MINV>
                    else
                        c1 = cmv(c2, mu);
                    end
                    c0 = dot(c1, mu, 1);
                end
                
                G.coef0 = c0;
                G.coef1 = c1;
                G.coef2 = c2;
                G.ldcov = ldc;
                G.has_cp = true;                
            end
        end        
        
        function G = from_cp(c1, c2, c0, ump)
            % Create Gaussian distribution(s) from canonical parameters
            %
            %   G = gaussd.from_cp(c1, c2);
            %   G = gaussd.from_cp(c1, c2, c0);
            %       creates an object to represent Gaussian distributions
            %       using canonical parameterization.
            %
            %       Input:
            %       - c1:   the coefficient of 1st order term [d x n matrix]
            %               (inv(sigma) * mu).
            %               For the models with zero means, c1 can be
            %               simply input as a scalar 0.
            %
            %       - c2:   the coefficient of 2nd order term
            %               (inv(sigma))
            %       - c0:   the constant term. The function will compute
            %               c0 if it is not specified.
            %
            %   G = gaussd.from_cp(c1, c2, [], 'mp');
            %   G = gaussd.from_cp(c1, c2, c0, 'mp');
            %       additionally derives the mean parameters.
            %
            
            % verify input types
            
            if ~(isfloat(c1) && ndims(c1) == 2)
                error('gaussd:from_cp:invalidarg', ...
                    'c1 should be a numeric matrix.');
            end
            
            if ~isobject(c2)
                error('gaussd:from_cp:invalidarg', ...
                    'c2 should be an object of symmetric matrix.');
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
            
            if isequal(c1, 0)
                d = c2.d;
                n = c2.n;             
            else
                [d, n] = size(c1);
                if c2.d ~= d
                    error('gaussd:from_cp:invalidarg', ...
                        'The dimension of c2 does not match that of c1.');
                end
                
                sn = c2.n;
                if sn > 1 && sn ~= n
                    error('gaussd:from_cp:invalidarg', ...
                        'The c2.n does not match size(c1, 2).');
                end
            end
            
            
            % determine mu and sigma                       
            
            if ump                
                sigma = inv(c2);
                if isequal(c1, 0)
                    mu = 0;
                else
                    if sn == 1
                        mu = sigma * c1; %#ok<MINV>
                    else
                        mu = cmv(sigma, c1);
                    end
                end                
            else
                if isequal(c1, 0)
                    mu = 0;
                else
                    if sn == 1
                        mu = c2 \ c1;
                    else
                        mu = cdv(c2, c1);
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
                if isequal(c1, 0)
                    c0 = 0;
                else
                    c0 = dot(c1, mu, 1);
                end
            end                      
            
            ldc = -lndet(c2);
            
            % construct Gaussian object
            
            G = gaussd();
            
            G.dim = d;
            G.num = n;
            G.coef0 = c0;
            G.coef1 = c1;
            G.coef2 = c2;
            G.ldcov = ldc;
            
            G.has_cp = true;
            
            if ump
                G.mean = mu;
                G.cov = sigma;
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
            %   D = sqmahdist(G, X);
            %       computes the squared Mahalanobis distances between the
            %       samples in X and the centers of the Gaussian 
            %       distributions selected by indices si.
            %
            %   Note that canonical parameters are required for the 
            %   computation.
            %
            
            if ~G.has_cp
                error('gaussd:sqmahdist:nocp', 'Canonical parameters are required.');
            end
            
            % take parameters
            
            c0 = G.coef0;
            c1 = G.coef1;
            c2 = G.coef2;
            
            if nargin >= 3
                c0 = c0(:, si);
                if ~isequal(c1, 0)
                    c1 = c1(:, si);
                end
                if c2.n ~= 1
                    c2 = c2.take(si);
                end
            end                        
                
            % compute 
            
            t2 = quad(c2, X, X);
            
            if isequal(c1, 0)
                D = t2;
            else
                t1 = c1' * X;
                
                n1 = size(c1, 2);
                n2 = c2.n;
                if n1 == 1
                    D = t2 - 2 * t1 + c0;
                else
                    if n2 == 1
                        D = bsxfun(@minus, t2, 2 * t1);
                    else
                        D = t2 - 2 * t1;
                    end
                    D = bsxfun(@plus, D, c0.');
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
                if G.coef2.n == 1
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
        
        
        function L = pdf(G, X, si)
            % Compute the probability density function
            %                        
            %   L = pdf(G, X)
            %       compute probability density function at the samples 
            %       given by columns of X.
            %
            %       Let m be the number of distributions contained in G,
            %       and n be the number of samples in X. Then L will be
            %       a matrix of size m x n, with L(i, j) being the pdf
            %       value at X(:,j) w.r.t the j-th Gaussian.
            %
            %   L = pdf(G, X, si);
            %       compute the probability density with respect to the 
            %       models selected by the index vector si.
            %
            %   Note that canonical parameters are required for the 
            %   computation.  
            %
            
            if nargin < 3
                L = exp(logpdf(G, X));
            else
                L = exp(logpdf(G, X, si));
            end            
        end
               
        
        function Gp = posterior(G, c1a, c2a, ump)
            % Get the posterior Gaussian models
            %
            %   Gp = posterior(G, c1a, c2a);
            %       compute the posterior Gaussian model(s) with G
            %       regarded as the prior.
            %
            %       c1a and c2a are quantities summarized from the
            %       observations, which are to be added to coef1 and
            %       coef2 respectively.                        
            %       
            %   Gp = posterior(G, c1a, c2a, 'mp');
            %       also derive the mean parameters
            %
            
            c10 = G.coef1;
            c20 = G.coef2;
            
            if size(c10, 2) == size(c1a, 2)
                c1 = c10 + c1a;
            else
                c1 = bsxfun(@plus, c10, c1a);
            end
            
            c2 = c20 + c2a;
            
            if nargin <= 3
                Gp = gaussd.from_cp(c1, c2);
            else
                if ~(ischar(ump) && strcmp(ump, 'mp'))
                    error('gaussd:posterior:invalidarg', ...
                        'the 4th argument can only be ''mp''.');
                end
                Gp = gaussd.from_cp(c1, c2, [], 'mp');
            end                        
        end
        
        
        function mu = pos_mean(G, c1a, c2a)
            % Compute the mean of the posterior Gaussian
            %
            %   mu = pos_mean(G, c1a, c2a);
            %       computes the mean of posterior Gaussian model
            %       with G regarded as the prior.
            %
            %       c1a and c2a are quantities summarized from the 
            %       observations, which are to be added to coef1 and
            %       coef2 respectively.
            %
            %   Note that G must contain only one distribution. 
            %
            
            if G.num ~= 1
                error('gaussd:pos_mean:invalidarg', ...
                    'The object must contain exactly one distribution.');
            end
            
            c1 = G.coef1 + c1a;
            c2 = G.coef2 + c2a;
            
            mu = c2 \ c1;
        end
                        
    end    
        
end



