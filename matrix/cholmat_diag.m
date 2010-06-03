classdef cholmat_diag < cholmat
    % The class of cholesky matrices by decomposing diagonal matrices
    %
    
    % Created by Dahua Lin, on Mar 22, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties
        s;      % the value of diagonal elements
    end
    
    methods
        function obj = cholmat_diag(s)
            % Constructs a matrix representing the result of Cholesky
            % decomposition of diagonal matrices.
            %
            %   obj = cholmat_diag(d, s);
            %       Constructs a matrix representing the result of
            %       Cholesky decomposition of diagonal matrices.
            %            
            %       s should be a matrix of size d x m, to 
            %       represent m matrices of size d x d.
            %
                       
            if ~(isfloat(s) && isreal(s) && ndims(s) == 2)
                error('cholmat_udiag:invalidarg', ...
                    's should be either a real scalar or a real row vector.');
            end
            
            obj = obj@cholmat(size(s, 1), size(s, 2));
            obj.s = s;
        end
        
        function o = take(obj, i)
            % Take a subset of the matrices contained in the object
        
            o = cholmat_diag(obj.s(:, i));        
        end
        
        
        function A = fullform(obj)
            % Gets the full form of size d x d x n.
                        
            A = makediag(obj.s);            
        end
        
        
        function Y = mtimes(obj, X)
            % Transforms the samples in X (compute C * X)
        
            if obj.num ~= 1
                error('cholmat_udiag:mtimes:nomultmat', ...
                    'mtimes only applies when obj.num == 1.');
            end
            
            if obj.dim == 1
                Y = obj.s * X;
            else
                Y = bsxfun(@times, obj.s, X);
            end
        end
                
    end
end