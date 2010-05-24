classdef pdmat
    % The base class for positive or semi-positive definite matrices
    %
    %   A positive (semi)definite matrix can be represented by
    %   a d x d array in MATLAB, this way is often not very
    %   efficient when the matrix assumes special structures.
    %
    %   Under specific conditions, more concise and efficient
    %   representation is available, such as when different
    %   elements are independent, etc. 
    %
    %   The pdmat base class defines some basic interfaces that
    %   a positive definite matrix should support.
    %   By implementing a derived class, one can create his
    %   own representation of positive matrix tailored to 
    %   specific needs, while maintaining the ability of
    %   interacting with other matrices and vectors.
    %
    %   This series of classes are designed to faciliate the
    %   development of models that rely on positive definite
    %   matrices (such as Gaussian models) such that they
    %   can be written in a generic way while enjoying the
    %   efficiency by delegating some core calculations to
    %   the derived class of pdmat.
    %   
    %   Each object can contain one or more actual matrices.
    %
    
    % Created by Dahua Lin, on Mar 22, 2010
    % Modified by Dahua Lin, on Apr 14, /afs/csail.mit.edu/u/d/dhlin/works/smitoolbox/gauss2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        dim;        % the dimension of the underlying vector space (d)
        num;        % the number of actual matrices contained in the object        
    end
    
    methods
        function obj = pdmat(d, m)
            % constructs an object
            %
            %   obj = pdmat(d, m);
            %       constructs an object with m matrices of size d x d.            
            %
            
            obj.dim = d;            
            obj.num = m;
        end
    end
    
    methods(Abstract)                
        R = take(obj, i);
        % takes a subset of the matrices contained in the object
        % i gives the index (indices) of them
        
        R = inv(obj);
        % compute the inverse of all matrices 
        
        b = is_subtype(obj, obj2);
        % tests whether obj2 is a special case of the class of this obj.
        
        Y = cmult(obj, X);
        % compute the product with corresponding vectors
        %
        % If there are m models, then X should be of size d x m.
        % In the output, Y is of size d x m, such that Y(:,i)
        % equals the multiplication of the i-th matrix and X(:,i).
        %
        
        Y = cldiv(obj, X)
        % compute the left matrix division with corresponding vectors
        %
        % If there are m models, then X should be of size d x m.
        % In the output, Y is of size d x m, such that Y(:,i)
        % equals M_i \ X(:,i), where M_i is the i-th matrix in obj.
        %
        
        y = mtimes(x1, x2);
        % compute the matrix multiplication with X
        %
        %   c * obj:    multiply a scalar c with the matrix
        %               returns an object of the same class
        %
        %               There are three different cases:
        %               - c is a scalar, and obj.num == m
        %                 multiplies c with each matrix
        %               - c is 1 x m, and obj.num == m
        %                 multiplies c(i) with the i-th matrix
        %               - c is 1 x m, and obj.num == 1
        %                 multiplies c(i) with the matrix, creating
        %                 different scaled version of the current matrix.
        %
        %               In all these cases, the output has 
        %               y.dim == obj.dim and y.num == m.
        %
        %   obj * X:    matrix product with a d x n matrix X
        %               returns a numeric matrix of size d x n.
        %               This applies only to single-matrix object.
        %        
        
        R = plus(obj, rhs);
        % compute the addition of two matrices
        %
        % The number of models in obj and rhs can be in either of the
        % following combination:
        %   obj.num == 1 && rhs.num == 1
        %   obj.num == 1 && rhs.num == m  (obj added to each matrix of rhs)
        %   obj.num == m && rhs.num == 1  (rhs added to each matrix of obj)
        %   obj.num == m && rhs.num == m  (addition done correspondingly)
        %  
        % Note that rhs should be of a sub-type of obj, or vice versa.
        %
        
        Y = mldivide(obj, X)
        % compute inv(M) * X
        %
        % Here, obj should have obj.num == 1 
        % and X is a matrix of size d x n, n can be any positive integer.
        %
        
        R = combine(obj, c);
        % Make a single-matrix object by pooling all matrices in the 
        % object
        % 
        % The output object is formed by c(1) * M_1 + ... + c(m) * M_m
        %
                                        
        Q = quadterm(obj, X);
        % compute the quadratic term on given samples
        %
        %   The sample matrix X should be of size d x n.
        %   If obj contains m matrices, then Q is a matrix of size m x n,
        %   where
        %
        %   Q(i, j) = X(:,j)' * C_i * X(:,j);
        %
        %   Here, C_i is the i-th matrix contained in the object.
        %
        
        R = fullform(obj);
        % convert the matrices to full d x d x m form
        
        R = chol(obj);
        % compute the (lower) Cholesky transform of all matrices
        %
        % It returns an object of class cholmat
        %
        
        R = sqrtm(obj);
        % compute the square root of all matrices
        %
                                
        v = logdet(obj);
        % compute the logarithm of determinants of all matrices
        %
        % In the output, v is of size 1 x m.    
        
        v = trace(obj);
        % compute the trace of all matrices
        %
        % In the output, v is of size 1 x m.
        %
    end
    
end

