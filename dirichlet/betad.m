classdef betad < expfdistr
    % The class to represent a beta distribution
    %
    %   The probability density function of a beta distribution
    %   parameterized by alpha and beta is 
    %
    %   p(x) = 1/B(alpha, beta) * x^(alpha-1) * (1-x)^(beta-1).
    %
    %   Each object can contain one or multiple distributions
    %
    
    % Created by Dahua Lin, on Mar 21, 2010
    % Modified by Dahua Lin, on April 14, 2010
    %
    
    properties
        alpha;      % alpha-parameter [1 x m]    
        beta;       % beta-parameter [1 x m]
    end
    
    properties(Dependent)
        mean;
        mode;
        var;
    end
    
    methods
        function v = get.mean(obj)
            a = obj.alpha;
            b = obj.beta;
            v = a ./ (a + b);
        end
        
        function v = get.mode(obj)
            a = obj.alpha;
            b = obj.beta;
            v = (a - 1) ./ (a + b - 2);
            v(a <= 1 & b <= 1) = nan;  % only applies to when a > 1 and b > 1
        end
        
        function v = get.var(obj)
            a = obj.alpha;
            b = obj.beta;
            s = a + b;
            v = (a .* b) ./ (s.^2 .* (s + 1)); 
        end
    end
    
    methods
        function obj = betad(alpha, beta)
            % constructs an object of beta distribution(s).
            %
            %   obj = betad(alpha, beta);
            %       constructs an object that contain beta distribution(s).
            %
            %       If there are m distributions, then both alpha and
            %       beta should be 1 x m vectors.
            %
            
            if ~(isfloat(alpha) && isreal(alpha) && ...
                    ndims(alpha) == 2 && size(alpha, 1) == 1)
                error('betad:invalidarg', 'alpha should be a real row vector.');
            end
            
            if ~(isfloat(beta) && isreal(beta) && ...
                    ndims(beta) == 2 && size(beta, 1) == 1)
                error('betad:invalidarg', 'beta should be a real row vector.');
            end
            
            if size(alpha, 2) ~= size(beta, 2)
                error('betad:invalidarg', ...
                    'The sizes of alpha and beta are not the same.');
            end
        
            obj.nmodels = size(alpha, 2);            
            obj.alpha = alpha;
            obj.beta = beta;
            
            obj.logpar = gammaln(alpha) + gammaln(beta) - gammaln(alpha + beta);            
        end
        
        function LT = compute_clinterm(obj, X, i)
            % Evaluate the canonical linear term on given samples
            %
            %   LT = compute_clinterm(obj, X)
            %   LT = compute_clinterm(obj, X, i);
            %       compute the canonical linear term on the 
            %       samples given in X for the models selected by
            %       index i.
            %       
            %       If i is omitted, then all models are selected.
            %
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == 1)
                error('betad:compute_clinterm:invalidarg', ...
                    'X should be a row vector.');
            end
            
            a = obj.alpha;
            b = obj.beta;
            if nargin >= 3 && ~isempty(i)
                a = a(i);
                b = b(i);
            end
            
            LT = (a - 1)' * log(X) + (b - 1)' * log(1 - X);
        end
        
        function L = compute_logbase(obj, X) %#ok<INUSD,MANU>
            % compute the logarithm of base measure on given samples
            L = 0;
        end
    end
    
end

