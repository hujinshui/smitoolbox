classdef gaussd_mp 
    % The class to represent Gaussian distribution with mean parameters
    %
    % In general, the mean parameterized form is not efficient for
    % evaluating pdf. One can use the to_cp method to convert a
    % mean-parameterized object to the corresponding object 
    % parameterized by canonical parameters.
    %
    
    % Created by Dahua Lin, on Mar 23, 2010
    %
    
    properties(GetAccess='public', SetAccess='protected')
        dim;        % the underlying space dimension (d)
        nmodels;    % the number of models (m)
        mu;         % the parameters corresponding to mean vectors [d x m]
        sigma;      % the parameters corresponding to variance/covariance (pdmat object)
    end
    
    methods
        
        function obj = gaussd_mp(mu, sigma)
            % Constructs a Gaussian distribution object with mean
            % parameters
            %
            %   obj = gaussd_mp(mu, sigma);
            %       constructs an object to represent Gaussian 
            %       distribution with mean parameterization.
            %
            %       If there are m distributions in a d-dimensional
            %       space, then mu should be a matrix of size d x m,
            %       and sigma should be a pdmat object with dim == d
            %       and num == m or num == 1 (when covariance is shared).            
            %   
            
            if ~(isfloat(mu) && ndims(mu) == 2)
                error('gaussd_mp:invalidarg', 'mu should be numeric matrix.');
            end
            
            [d, m] = size(mu);
            
            if ~(isa(sigma, 'pdmat') && sigma.dim == d && ...
                (sigma.num == 1 || sigma.num == m))
                error('gaussd_mp:invalidarg', ...
                    'sigma should be a pdmat object with dim == d and num being 1 or m.');
            end
            
            obj.dim = d;
            obj.nmodels = m;
            obj.mu = mu;
            obj.sigma = sigma;
        end
                
        function R = to_cp(obj)
            % Get the corresponding canonical parameterized object
            %
            %   R = to_cp(obj);
            %       returns the Gaussian distribution object with
            %       canonical parameterization
            %
            
            R = gaussd_cp.from_mean_and_cov(obj.mu, obj.sigma);
        end
        
        function S = get_sampler(obj, rstream)
            % Get the sampler for drawing samples from the distribution
            %
            %   S = obj.get_sampler();
            %   S = obj.get_sampler(rstream);
            %       Gets the sampler S for drawing samples from the
            %       distribution
            %
            %       rstream is a random strea object as random source.
            %       If it is not specified, the default stream is used.
            %
            
            if nargin < 2
                rstream = RandStream.getDefaultStream();
            end
            
            S = gaussd_sampler(obj, rstream);                
        end
        
    end
end
