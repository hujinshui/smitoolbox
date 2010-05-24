classdef cholmat_gen < cholmat
    % The class of cholesky matrices by decomposing generic matrices    
    %
    
    % Created by Dahua Lin, on Mar 22, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties
        mats;      % the matrices
    end
    
    methods
        function obj = cholmat_gen(C)
            % Constructs a matrix representing the result of Cholesky
            % decomposition of a uniform diagonal matrix.
            %
            %   obj = cholmat_gen(C);
            %       Constructs a matrix representing the result of
            %       Cholesky decomposition of generic matrices.            
            %
            %       C is an array of size d x d x m, where each C(:,:,i)
            %       is a lower triangle matrix.
            %
           
            if ~(isfloat(C) && isreal(C) && ndims(C) <= 3 && size(C,1) == size(C,2))
                error('cholmat_gen:invalidarg', ...
                    'C should be a d x d x m numeric array.');
            end
                                    
            obj = obj@cholmat(size(C,1), size(C,3));
            obj.mats = C;
        end
        
        function o = take(obj, i)
            % Take a subset of the matrices contained in the object
        
            o = cholmat_gen(obj.mats(:,:,i));        
        end
        
        
        function A = fullform(obj)
            % Gets the full form of size d x d x n.
            
            A = obj.mats;
        end
        
        
        function Y = mtimes(obj, X)
            % Transforms the samples in X (compute C * X)
        
            if obj.num ~= 1
                error('cholmat_gen:mtimes:nomultmat', ...
                    'mtimes only applies when obj.num == 1.');
            end
            
            Y = obj.mats * X;
        end
                
    end
end