classdef gaussgm
    % The class for inference/estimation on Gaussian generative model
    %
    % Remarks
    % -------
    %   - The parameters of this model is a gaussd object.
    %   - This class adapts Gaussian distribution to a generative model
    %     without prior. This is useful in cooperating with fmm.
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Nov 14, 2010
    %
    
    
    properties(GetAccess='public', SetAccess='private')
        dim;
    end
    
    properties        
        cov_type = 'full';
        tie_cov = false;        
    end
    
    
    methods
        
        function obj = set.cov_type(obj, v)
            if ~(ischar(v))
                error('gaussgm:invalidprop', ...
                    'cov_type should be a string.');
            end
            if strcmp(v, 'full') || strcmp(v, 'diag') || strcmp(v, 'iso')
                obj.cov_type = v;
            else
                error('gaussgm:invalidprop', ...
                    'cov_type can be either of ''full'', ''diag'', or ''iso''.');
            end
        end
        
        
        function obj = set.tie_cov(obj, v)            
            if ~(islogical(v) && isscalar(v))
                error('gaussgm:invalidprop', ...
                    'tie_cov should be a logical scalar.');
            end
            obj.tie_cov = v;            
        end
        
    end
    
    
    methods
        
        function obj = gaussgm(d)
            % Constructs a Gaussian generative model
            %
            %   obj = gaussgm(d);
            %       constructs a Gaussian generative model on a
            %       vector space of dimension d.
            %
            
            if ~(isscalar(d) && isnumeric(d) && d == fix(d) && d >= 1)
                error('gaussgm:invalidarg', 'd should be a positive integer.');
            end            
            obj.dim = d;
        end
        
        
        function n = check_parameters(obj, g)
            % Check validity of parameters and return the number
            %
            %   n = obj.check_parameters(g);
            %       if g is a valid gaussd object, it returns
            %       the number of models, otherwise it returns -1.
            %
            
            if isa(g, 'gaussd') && g.dim == obj.dim
                n = g.num;
            else
                n = -1;
            end
            
        end
        
        
        function n = check_observations(obj, X)
            % Check validity of observations and return the number
            %
            %   n = obj.check_observations(X);
            %       if X is a valid observation matrix, it returns
            %       the number of samples in X, otherwise it returns -1.
            %
            
            if isfloat(X) && ndims(X) == 2 && size(X,1) == obj.dim
                n = size(X,2);
            else
                n = -1;
            end            
        end                        
        
    end
    
                  
    
    methods                
        
        function LP = logpri(obj, g) %#ok<MANU>
            % Compute the log-prior of parameters
            %
            %   LP = obj.logpri(g);
            %
            
            LP = zeros(1, g.num);
        end
                                  
        
        function LL = loglik(obj, g, X) %#ok<MANU>
            % Compute log-likelihood
            %
            %   LL = obj.loglik(g, X);
            %       computes the log-likelihood at observed samples
            %       given by X with respect to the models g.
            %
            %       Suppose there are m models, and n samples, then
            %       LL will be a matrix of size m x n.
            %
            
            LL = g.logpdf(X);
        end                                                    
        
        function g = estimate_map(obj, X, w)
            % Performs MAP estimation
            %
            %   g = obj.estimate_map(X, w);
            %       
            %   Since no prior is used, MAP estimation is equivalent
            %   to MLE estimation.
            %
                        
            if nargin < 3
                w = 1;
            end
            
            g = gaussmle(X, w, ...
                'cov_type', obj.cov_type, ...
                'tie_cov', obj.tie_cov, ...
                'use_ip', true);
        end                        
                               
    end
    
end



