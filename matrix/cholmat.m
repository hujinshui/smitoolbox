classdef cholmat
    % The class for representing the matrices obtained from Cholesky
    % decomposition
    %
    
    % Created by Dahua Lin, on Mar 22, 2010
    %
    
    properties
        dim;        % the space dimension
        num;        % the number of matrices contained in the object
    end
    
    methods
        
        function obj = cholmat(d, m)
            % Constructing a matrix representing the Cholesky result
            %
            %   obj = cholmat(d, m);
            %       In the input, 
            %       d: the dimension, each matrix is d x d.
            %       m: the number of matrices contained in the object
            %
            
            obj.dim = d;
            obj.num = m;
        end
        
    end
    
    methods(Abstract)
        
        o = take(obj, i);
        % Take a subset of the matrices contained in the object
        %
        
        A = fullform(obj);
        % Gets the full form of size d x d x n.
        %
        
        Y = mtimes(obj, X);
        % Transforms the samples in X (compute C * X)
        %
        % Only works when obj.num == 1
        %
        
    end
end