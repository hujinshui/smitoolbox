classdef gausspri < prior_base
    % Gaussian prior distribution
    %
    %   This class wrap a gaussd object to prior_base interface
    %
    
    % Created by Dahua Lin, on Dec 27, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')
        dim;        % the space dimension
        gdistr;     % the Gaussian distribution (gaussd struct) 
        
        const_a;    % the Gaussian constant (a)
        const_b;    % the Gaussian constant (b)
    end
    
    methods
        
        function n = query_samples(obj, X)
            % Verify the validity of input samples and return the number
            %
            %   n = obj.query_samples(X);
            %       verifies the validity of X as a sample matrix, and
            %       returns the number of samples in X.
            %
            
            d = obj.dim;            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == d)
                error('query_samples:invalidarg', ...
                    'X should be a real matrix with size(X,1) == d.');
            end            
            n = size(X, 2);
        end
        
        
        function L = logpdf(obj, X)
            % Evaluate the log pdf at given samples
            %
            %   L = obj.logpdf(X);
            %       evaluates the log pdf at the samples in X.            
            %
            
            g = obj.gdistr;
            ca = obj.const_a;
            cb = obj.const_b;
            
            L = gaussd_logpdf(g, X, {ca, cb});            
        end
        
        
        function X = sample(obj, n)
            % Samples from the prior distribution
            %
            %   X = obj.sample(n);
            %       draws n samples from the Gaussian prior
            %
            
            g = obj.gdistr;
            X = gaussd_sample(g, n);
        end
        
        function X = pos_sample(obj, S, n)
            % Samples from posterior distribution 
            %
            %   X = obj.pos_sample(S, n);
            %       draws n samples from the posterior Gaussian
            %       distribution. 
            %
            %       Here, S is a gaussd struct that captures the 
            %       observation statistics. Note that S should have
            %       S.n == 1.
            %
            
            if S.n ~= 1
                error('gausspri:invalidarg', ...
                    'S violates the constraint: S.n == 1.');
            end            
            gp = gaussd_conjupdate(obj.gdistr, S);
            X = gaussd_sample(gp, n);
        end
        
        function X = mapest(obj, S)
            % Performs MAP estimation with the stats of observations
            %
            %   X = obj.mapest(S);
            %       performs MAP estimation with S capturing the 
            %       statistics from the observations.
            %
            
            g = obj.gdistr;
            X = gaussd_mapest(g, S);
        end
        
    end
    
end


