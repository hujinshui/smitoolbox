classdef cholmat_udiag < cholmat
    % The class of cholesky matrices by decomposing uniform diagonal
    % matrices
    %
    
    % Created by Dahua Lin, on Mar 22, 2010
    %
    
    properties
        s;      % the value of diagonal elements
    end
    
    methods
        function obj = cholmat_udiag(d, s)
            % Constructs a matrix representing the result of Cholesky
            % decomposition of a uniform diagonal matrix.
            %
            %   obj = cholmat_udiag(d, s);
            %       Constructs a matrix representing the result of
            %       Cholesky decomposition of a uniform diagonal matrix.
            %
            %       d is the dimension of the space, and s is the 
            %       value of the diagonal elements.
            %
            %       s should be a scalar or a 1 x m 
           
            if ~(isscalar(d) && d >= 1)
                error('cholmat_udiag:invalidarg', ...
                    'd should be a positive integer scalar.');
            end
            
            if ~(isfloat(s) && isreal(s) && ndims(s) == 2 && size(s, 1) == 1)
                error('cholmat_udiag:invalidarg', ...
                    's should be either a real scalar or a real row vector.');
            end
            
            obj = obj@cholmat(d, size(s, 2));
            obj.s = s;
        end
        
        function o = take(obj, i)
            % Take a subset of the matrices contained in the object
        
            o = cholmat_udiag(obj.dim, obj.s(i));        
        end
        
        
        function A = fullform(obj)
            % Gets the full form of size d x d x n.
            
            d = obj.dim;
            if d == 1
                A = makediag(obj.s);
            else
                A = makediag(obj.s(ones(d, 1), :));
            end        
        end
        
        
        function Y = mtimes(obj, X)
            % Transforms the samples in X (compute C * X)
        
            if obj.num ~= 1
                error('cholmat_udiag:mtimes:nomultmat', ...
                    'mtimes only applies when obj.num == 1.');
            end
            
            Y = obj.s * X;
        end
                
    end
end