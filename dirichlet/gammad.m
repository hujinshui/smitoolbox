classdef gammad < expfdistr
    % The class to represent Gamma distribution 
    %
    % A gamma distribution is parameterized by 
    %   - alpha:    the shape parameter
    %   - lambda:   the rate parameter (1 / scale parameter)
    % The probability density function is given by
    %
    %   p(x) = x^(alpha-1) * lambda^alpha * exp(-lambda * x) / Gamma(alpha)   
    %
    % The rate parameter lambda can be shared among all models contained
    % in the same object of gammad.
    %
    
    % Created by Dahua Lin, on Mar 21, 2010
    % Modified by Dahua Lin, on April 13, 2010
    %
    
    %% parameters
    
    properties(GetAccess='public', SetAccess='private')
        alpha;      % the shape parameter [1 x m]
        lambda;     % the rate parameter (1 / scale parameter) [1 x m or scalar]
    end
        
    %%  statistics of the distribution
    
    properties(Dependent)
        mean;
        mode;
        var;
    end
    
    methods
        function v = get.mean(obj)
            v = obj.alpha .* (1 ./ obj.lambda);
        end
        
        function v = get.mode(obj)
            v = (obj.alpha - 1) .* (1 ./ obj.lambda);
        end
        
        function v = get.var(obj)
            v = obj.alpha .* ((1 ./ obj.lambda) .^ 2);
        end
    end    
    
    
    %% main methods
    
    methods
        
        function obj = gammad(alpha, lambda)
            % Constructs an object for gamma distributions
            %
            %   obj = gammad(alpha, lambda);
            %       constructs an object of class gammad which contains
            %       one or multiple gamma distribution(s).
            %
            %       If there are m distributions, then 
            %       the shape parameter alpha should be a 1 x m vector,
            %       and the rate parameter lambda can be either a
            %       scalar (if all models shared the same rate) or
            %       a 1 x m vector otherwise.
            %
            
            if ~(isfloat(alpha) && ndims(alpha) == 2 && size(alpha,1) == 1)
                error('gammad:invalidarg', ...
                    'alpha should be a 1 x m numeric vector.');
            end
            
            m = size(alpha, 2);
            
            if ~(isfloat(lambda) && ndims(lambda) == 2 && size(lambda, 1) == 1 && ...
                (size(lambda, 2) == 1 || size(lambda, 2) == m))
                error('gammad:invalidarg', ...
                    'lambda should be either a scalar or of size 1 x m.');
            end
            
            obj.nmodels = m;
            
            obj.alpha = alpha;
            obj.lambda = lambda;
            
            obj.logpar = gammaln(alpha) - alpha .* log(lambda);
        end
        
        
        function LT = compute_clinterm(obj, X, i)
            % Evaluate the canonical linear term of the selected models on 
            % given samples          
            %
            
           if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == 1)
                error('gammad:compute_clinterm:invalidarg', ...
                    'X should be a numeric row vector.');
           end
            
            a = obj.alpha;
            r = obj.lambda;
            
            if nargin >= 3 && ~isempty(i)
                a = a(i);
                if ~isscalar(r)
                    r = r(i);
                end
            end                
            
            if size(a, 2) == size(r, 2)
                LT = (a - 1)' * log(X) - r' * X;
            else
                LT = bsxfun(@minus, (a-1)' * log(X), r' * X);
            end
        end
           
        function L = compute_logbase(obj, X) %#ok<INUSD,MANU>
            % compute the logarithm of base measure on given samples            
            L = 0;
        end
    end
    
end

