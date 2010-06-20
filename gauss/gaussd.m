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
                    'mu should be a two-dimensional matrix.');
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
            
        end
        
        
        
    end
    
    
    
end