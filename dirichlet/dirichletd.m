classdef dirichletd < expfdistr
    % The class to represent a Dirichlet distribution
    %
    %   A Dirichlet distribution is a distribution defined on 
    %   a (n-1)-dimensional probability simplex. It is parameterized
    %   by a n-dimensional vector alpha (concentration parameter).
    %   
    %   The probability density function is defined by 
    %
    %   p(x) = 1 / B(alpha) * \prod_{i=1}^d x_i^{alpha_i - 1}
    %
    %   Here, x_1 + ... + x_d = 1.
    %
    %   Each object of this class can contain one or multiple 
    %   dirichlet distributions defined of same dimensionality.
    %
    
    % Created by Dahua Lin, on Mar 21, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        dim;        % the dimension of the space (d)
        alpha;      % the concentration parameters [d x m]
        
        alpha_t;    % the sum of alpha values (total concentration) [1 x m]
    end
    
    properties(Dependent)
        mean;       % the mean value of x for each model [d x m]
    end    
    
    methods
        function v = get.mean(obj)
            v = bsxfun(@times, obj.alpha, 1 ./ obj.alpha_t);
        end
    end
            
    methods        
        function obj = dirichletd(alpha)
            % constructs a object of Dirichlet distributions
            %
            %   obj = dirichletd(alpha);
            %       constructs an object containing Dirichlet
            %       distributions.
            %
            %       If there are m models defined on a d-dimensional 
            %       space ((d-1)-simplex), then alpha should be a
            %       matrix of size d x m.
            %
            
            if ~(isfloat(alpha) && ndims(alpha) == 2)
                error('dirichletd:invalidarg', ...
                    'alpha should be a numeric matrix.');
            end
            
            [d, m] = size(alpha);
            
            obj.dim = d;
            obj.nmodels = m;
            
            % compute log-partition
            
            s = sum(alpha, 1);
            obj.alpha = alpha;
            obj.alpha_t = s;
            obj.logpar = sum(gammaln(alpha), 1) - gammaln(s);            
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
            %   The caller should ensure that each column of X 
            %   sums to unity.
            %
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == obj.dim)
                error('dirichletd:invalidarg', ...
                    'X should be a numeric matrix of size d x n.');
            end
            
            a = obj.alpha;
            if nargin >= 3 && ~isempty(i)
                a = a(:, i);
            end
            
            LT = (a - 1)' * log(X);                        
        end
        
        function L = compute_logbase(obj, X) %#ok<INUSD,MANU>
            % compute the logarithm of base measure on given samples
            L = 0;
        end
        
    end
end
