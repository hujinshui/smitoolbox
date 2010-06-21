classdef gaussd
    % The class to represent Gaussian distribution
    %
    %   Each object of this class can contain one or multiple Gaussian
    %   distributions. A Gaussian distribution can be parameterized by
    %   either canonical parameterization or mean parameteration.
    %
    %   Canonical parameterization is given by 
    %   - c1:   the coefficient vector of first order term
    %   - c2:   the coefficient matrix of second order term
    %   - c0:   the constant term
    %
    %   Mean parameterization is given by
    %   - mu:       the mean vector 
    %   - sigma:    the covariance matrix
    %  
    %   These two types of parameterization are related to each other
    %   as follows:
    %   - c2 = inv(sigma)
    %   - c1 = c2 * mu
    %   - c0 = mu' * inv(sigma) * mu + lndet(sigma) + d * log(2*pi)
    %   Or, equivalently,
    %   - sigma = inv(c2)
    %   - mu = sigma * c1
    %
    %   For Gaussian distributions with zero mean, one can simply set
    %   c1 or mu to a scalar 0 despite the actual dimension of the space.
    %   
    %   Generally, the evaluation of Mahalanobis distance or probability
    %   density function, and the posterior updating relies on the
    %   canonical parameterization, while maximum likelihood estimation
    %   would lead to mean parameterization.
    %
    %   The caller can specify which parameterization to use (or both)
    %   in constructing the object.
    %
    
    % Created by Dahua Lin, on June 19, 2010
    %
    
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')        
        dim;    % the dimension of the vector space
        num;    % the number of models contained in the object
        
        use_cp; % whether canonical parameterization is used
        use_mp; % whether mean parameterization is used
        
        mu;     % the mean vector(s)
        sigma;  % the covariance matrix(ces)
        
        c0;     % the constant term in canonical parameterization
        c1;     % the coefficient of the 1st order term (inv(sigma) * mu)
        c2;     % the coefficient of the 2nd order term (inv(sigma)                                
    end
    
    
    %% constructor
    
    methods
        
        function G = gaussd(d, n, mu, sigma, c0, c1, c2)
            % Construct a Gaussian distribution object
            %
            %   G = gaussd(d, n, mu, sigma, c0, c1, c2);
            %       
            %   This function is not recommended for direct use. 
            %   One can invoke static functions to create gaussd
            %   objects, such as gaussd.from_mp or gaussd.from_cp.
            %   
            
            G.dim = d;
            G.num = n;
            
            G.use_cp = ~(isempty(c0) || isempty(c1) || isempty(c2));
            G.c0 = c0;
            G.c1 = c1;
            G.c2 = c2;
            
            G.use_mp = ~(isempty(mu) || isempty(sigma));
            G.mu = mu;
            G.sigma = sigma;                        
        end

    end
    
    
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
            %               specified as a scalar 0.
            %       - sigma: the covariance matrix. It should be an
            %                object of a symmetric matrix class.
            %                (udmat, dmat, gsymat, etc).
            %
            %   G = gaussd.from_mp(mu, sigma, 'cp');
            %       additionally derives the canonical parameterization.
            %                                    
            
            if ~(isfloat(mu) && ndims(mu) == 2)
                error('gaussd:from_mp:invalidarg', ...
                    'mu should be a numeric matrix.');
            end
            
            if ~(isobject(sigma))
                error('gaussd:from_mp:invalidarg', ...
                    'sigma should be an object of symmetric matrix.');
            end
            
            if isequal(mu, 0)
                d = sigma.d;
                n = sigma.n;
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
            
            if ~ucp
                G = gaussd(d, n, mu, sigma, [], [], []);
            else
                c2_ = inv(sigma);
                
                if isequal(mu, 0)
                    c1_ = 0;
                    c0_ = lndet(sigma) + d * log(2*pi);
                else
                    if sn == 1
                        c1_ = c2_ * mu; %#ok<MINV>
                    else
                        c1_ = cmv(c2_, mu);
                    end
                    c0_ = sum(c1_ .* mu, 1) + lndet(sigma) + d * log(2*pi);
                end
                G = gaussd(d, n, mu, sigma, c0_, c1_, c2_);
            end
        end
        
        
        function G = from_cp(c1, c2, c0)
            % Create Gaussian distribution(s) from canonical
            % parameterization
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
            %   Note that the mean parameterization will also be 
            %   derived during the construction.
            %
            
            if ~(isfloat(c1) && ndims(c1) == 2)
                error('gaussd:from_cp:invalidarg', ...
                    'c1 should be a numeric matrix.');
            end
            
            if ~isobject(c2)
                error('gaussd:from_cp:invalidarg', ...
                    'c2 should be an object of symmetric matrix.');
            end
            
            if isequal(c1, 0)
                d = c2.d;
                n = c2.n;
                sigma_ = inv(c2);
                mu_ = 0;                
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
                
                sigma_ = inv(c2);
                if sn == 1
                    mu_ = sigma_ * c1; %#ok<MINV>
                else
                    mu_ = cmv(sigma_, c1);
                end
            end
                                                            
            if nargin >= 3 && ~isempty(c0)
                if ~(isfloat(c0) && isequal(size(c0), [1 n]))
                    error('gaussd:from_cp:invalidarg', ...
                        'c0 should be a numeric vector of size 1 x n.');
                end
            else
                c0 = lndet(sigma_) + d * log(2*pi);
                if ~isequal(c1, 0)
                    c0 = sum(c1 .* mu_, 1) + c0;
                end
            end                      
                                
            G = gaussd(d, n, mu_, sigma_, c0, c1, c2);                
        end
    end
        
    
    %% Probability evaluation
    
    methods
        
        function L = logprob(G, X, si)
            % Compute logarithm of PDF of given samples
            %
            %   L = logprob(G, X)
            %       compute logarithm of probability density function
            %       at the samples given by columns of X.
            %
            %       Let m be the number of distributions contained in G,
            %       and n be the number of samples in X. Then L will be
            %       a matrix of size m x n, with L(i, j) being the pdf
            %       value at X(:,j) w.r.t the j-th Gaussian.
            %
            %   L = logprob(G, X, si);
            %       compute the logarithm of pdf with respect to the 
            %       models selected by the index vector si.
            %
            
            if ~G.use_cp
                error('gaussd:logprob:nocp', ...
                    'Canonical parameterization is required.');
            end
            
            c0_ = G.c0;
            c1_ = G.c1;
            c2_ = G.c2;
            n_ = G.num;
            
            if n_ == 1   % single model                
                if isequal(c1_, 0) 
                    t2 = quad(c2_, X, X);
                    L = -0.5 * (t2 + c0_);
                else
                    t1 = (-2) * c1_' * X;
                    t2 = quad(c2_, X, X);                    
                    L = -0.5 * (t1 + t2  + c0_);
                end                
                if nargin >= 3 && ~isequal(si, 1)
                    L = L(si, :);
                end                
                
            elseif c2_.n == 1 % tied covariance
                if nargin < 3
                    t1 = (-2) * c1_' * X;
                    t0 = c0_;
                else
                    t1 = (-2) * c1_(:, si)' * X;
                    t0 = c0_(si);
                end
                t2 = quad(c2_, X, X);
                if size(t1, 1) == 1
                    L = -0.5 * (t1 + t2 + t0);
                else
                    L = -0.5 * bsxfun(@plus, ...
                        bsxfun(@plus, t1, t2), t0.');
                end       
                
            else  % multiple models with respective covariance
                if isequal(c1_, 0)
                    if nargin < 3
                        t2 = quad(c2_, X, X);
                        t0 = c0_;
                    else
                        t2 = quad(c2_.take(si), X, X);
                        t0 = c0_(si);
                    end
                    L = -0.5 * bsxfun(@plus, t2, t0.');
                else
                    if nargin < 3
                        t1 = (-2) * c1_' * X;
                        t2 = quad(c2_, X, X);
                        t0 = c0_;
                    else
                        t1 = (-2) * c1_(:, si)' * X;
                        t2 = quad(c2_.take(si), X, X);
                        t0 = c0_(si);
                    end
                    L = -0.5 * bsxfun(@plus, t1 + t2, t0.');
                end
            end
            
        end
        
        
        function Gpos = inject(Gpri, dc1, dc2)
            % Get posterior Gaussian by incorporating updates to c1 and c2
            %
            %   Gpos = Gpri.inject(dc1, dc2);
            %       it returns a Gaussian distribution whose canonical
            %       parameters are given by c1 + dc1 and c2 + dc2.
            %
            
            Gpos = gaussd.from_cp(Gpri.c1 + dc1, Gpri.c2 + dc2);
        end
                        
    end
    
        
end



