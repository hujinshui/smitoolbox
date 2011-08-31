classdef gammad
    % The class to represent (one or multi-dimensional) Gamma distribution
    %
    %   Each object of this class can contain one or multiple Gamma
    %   distributions. A Gamma distribution is characterized by two
    %   parameters:
    %
    %   - shape:    the shape parameter 
    %   - scale:    the scale parameter
    %  
    %   For an object with n distributions over d-dimensional space:    
    %   - shape is a matrix of size 1 x n or d x n.
    %   - scale can be a 1 x n row vector or a scale (shared)    
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        
        dim;    % the dimension of the underlying space
        num;    % the number of distributions contained in the object
        
        shape;  % the shape parameter(s)
        scale;  % the scale parameter(s)
    end
        
    methods
        
        function v = mean(obj)
            % Evaluates the mean(s) of the distribution
            %
            %   v = mean(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding 
            %       to the i-th distribution in obj.
            %
        
            k = obj.shape;
            sc = obj.scale;
            d = obj.dim;
            
            if isscalar(sc)
                if sc == 1
                    v = k;
                else
                    v = k * sc;
                end
            else
                v = bsxfun(@times, k, sc);
            end
            
            if size(k, 1) ~= d
                v = v(ones(d, 1), :);
            end                
        end
                
        function v = var(obj)
            % Evaluates the variance(s) of the distribution
            %
            %   v = var(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding 
            %       to the i-th distribution in obj.
            %
            
            k = obj.shape;
            sc = obj.scale;
            d = obj.dim;
            
            if isscalar(sc)
                if sc == 1
                    v = obj.shape;
                else
                    v = (sc^2) * obj.shape;
                end
            else
                v = bsxfun(@times, obj.shape, (sc.^2));
            end
            
            if size(k, 1) ~= d
                v = v(ones(d, 1), :);
            end
        end
            
        
        function v = entropy(obj)
            % Evaluates the entropy value(s) of the distribution
            %
            %   v = entropy(obj);
            %
            %       v will be a 1 x n matrix, with v(:,i) corresponding 
            %       to the i-th distribution in obj.
            %            
            
            k = obj.shape;
            sc = obj.scale;
            d = obj.dim;
            
            v = k + gammaln(k) + (1 - k) .* psi(k);
            if ~isequal(sc, 1)
                if isscalar(sc) 
                    v = v + log(sc);
                else
                    v = bsxfun(@plus, v, log(sc));
                end
            end
            
            if size(k, 1) == 1
                if d > 1
                    v = v * d;
                end
            else
                v = sum(v, 1);
            end
        end
        
    end
    
    %% Construction
    
    methods 
        
        function obj = gammad(shape, scale, d)
            % Constructs a Gamma distribution object
            %
            %   obj = gammad(shape);
            %   obj = gammad(shape, scale);
            %
            %       constructs a Gamma distribution object given 
            %       the parameters.
            %
            %       Inputs:
            %       - shape:    the shape parameter of size d x n.
            %
            %       - scale:    the scale parameter, which can be
            %                   either a 1 x n row vector or a 
            %                   scale (if shared by all distributions).
            %
            %   obj = gammad(shape, scale, d);
            %   
            %       To construct multi-dimensional gamma distributions,
            %       where each dimension has the same shape param, then
            %       one can use this syntax.
            %
            %       Here, shape can be input a 1 x n row vector, and
            %       use the 3rd argument to specify the dimension.
            %
                        
            % verify inputs
            
            if ~(isfloat(shape) && ndims(shape) == 2)
                error('gammad:invalidarg', ...
                    'shape should be a numeric matrix.');
            end            
            [d_, n] = size(shape);
            
            if ~( isfloat(scale) && (isscalar(scale) || ...
                    isequal(size(scale), [1 n])) )
                error('gammad:invalidarg', ...
                    'scale should be a scalar or a 1 x n row vector.');
            end
            
            if nargin < 3
                d = d_;
            else
                if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1)
                    error('gammad:invalidarg', ...
                        'd should be a positive integer scalar.');
                end
                if d_ > 1 && d ~= d_
                    error('gammad:invalidarg', ...
                        'The size of shape is inconsistent with d.');
                end
            end
            
            obj.dim = d;
            obj.num = n;
            obj.shape = shape;
            obj.scale = scale;
        end            
    end
    
    
    %% Evaluation
    
    
    %% Sampling
    
    
    
end




