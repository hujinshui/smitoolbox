classdef pdmat_udiag < pdmat
    % The object to represent uniform diagonal Gaussian matrices with the 
    % same value at each diagonal elements (dv)
    %
    % It equals dv * eye(d).
    %
    
    % Created by Dahua Lin, on Mar 23, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties
        dv;       % the value of diagonal entries
    end
    
    methods
        function obj = pdmat_udiag(d, dv)
            % construct an object from the given diagonal value
            %
            %   obj = pdmat_udiag(d, dv);
            %       constructs an object that contains the d x d Gaussian
            %       matrices with uniform diagonal elements given by dv.
            %
            %       dv should be a 1 x m row vector, then the i-th            
            %       matrix is dv(i) * eye(d).
            %
                     
            if ~(isnumeric(d) && isscalar(d) && d >= 1)
                error('pdmat_udiag:invalidarg', 'd should be an integer scalar.');
            end
            
            if ~(isfloat(dv) && ndims(dv) == 2 && size(dv,1) == 1)
                error('pdmat_udiag:invalidarg', ...
                    'dv should be a scalar or a numeric row vector.');
            end
                                                
            obj = obj@pdmat(d, size(dv, 2));            
            obj.dv = dv;
        end
        
        
        function R = take(obj, i)
            % takes a subset of the matrices contained in the object
            
            R = pdmat_udiag(obj.dim, obj.dv(i));
        end

        
        function R = inv(obj)
            % compute the inverse of all matrices 
                                                
            R = pdmat_udiag(obj.dim, 1 ./ obj.dv);
        end
        
        
        function b = is_subtype(obj, obj2) %#ok<MANU>
            % tests whether obj2 is of a sub-type of the class pdmat_udiag.
            
            b = isa(obj2, 'pdmat_udiag');
        end
        
        
        function Y = cmult(obj, X)
            % compute the product with corresponding vectors
            
            d = obj.dim;
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [d m]))
                error('pdmat_udiag:cmult:invalidarg', ...
                    'X should be a d x m matrix.');
            end
                                    
            if d == 1
                Y = obj.dv .* X;
            else
                Y = bsxfun(@times, obj.dv, X);
            end
        end        
        
        
        function Y = cldiv(obj, X)
            % compute the left matrix division with corresponding vectors
            
            d = obj.dim;
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [d m]))
                error('pdmat_udiag:cldiv:invalidarg', ...
                    'X should be a d x m matrix.');
            end
                                    
            if d == 1
                Y = X ./ obj.dv;
            else
                Y = bsxfun(@times, (1 ./ obj.dv), X);
            end
        end
        
        
        function y = mtimes(x1, x2)
            % compute the multiplication
            
            if isa(x1, 'pdmat_udiag')
                if x1.num ~= 1
                    error('pdmat_udiag:mtimes:nomultmat', ...
                        'Multiplication with a multi-matrix object is not allowed.');
                end
            
                if ~(isfloat(x2) && ndims(x2) == 2 && size(x2,1) == x1.dim)
                    error('pdmat_udiag:mtimes:invalidarg', ...
                        'The right multiplicand should be a d x n matrix.');
                end
            
                y = x1.dv * x2;
            else
                if ~isfloat(x1) || ... 
                        ~(isscalar(x1) || x2.num == 1 || isequal(size(x1), [1 x2.num]))
                    error('pdmat_udiag:mtimes:invalidarg', ...
                        'The left multiplicand should be either a scalar or a 1 x m vector.');
                end
            
                y = pdmat_udiag(x2.dim, x1 .* x2.dv);
            end
        end
        
        
        function R = plus(obj, rhs)
            % compute the addition of two matrices
            
            if obj.dim ~= rhs.dim
                error('pdmat_udiag:plus:invalidarg', ...
                    'The dimensions of addends are not the same.');
            end
                                                            
            if is_subtype(obj, rhs)
                if ~(obj.num == 1 || rhs.num == 1 || obj.num == rhs.num)
                    error('pdmat_udiag:plus:invalidarg', ...
                        'Invalid number of matrices in the objects.');
                end
                
                R = pdmat_udiag(obj.dim, obj.dv + rhs.dv);
                
            elseif is_subtype(rhs, obj)
                
                R = plus(rhs, obj);
                
            else
                error('pdmat_udiag:incompatible', ...
                    'The addend is incompatible with the object.');
            end            
        end
        
        
        function Y = mldivide(obj, X)
            % Compute inv(M) * X
            
            if obj.num ~= 1
                error('pdmat_udiag:mldivide:momultmat', ...
                    'Do mldivide with a multi-matrix object is not allowed.');
            end
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == obj.dim)
                error('pdmat_udiag:mldivide:invalidarg', ...
                    'The right multiplicand should be a d x n numeric matrix.');
            end
            
            Y = X * (1 / obj.dv);            
        end
        
        
        function R = combine(obj, c)
            % Make a single-matrix object by combining all matrices
            
            m = obj.num;            
            if ~(isfloat(c) && isequal(size(c), [1 m]))
                error('pdmat_udiag:combine:invalidarg', ...
                    'c should be a 1 x m row vector.');
            end
            
            if m == 1
                if c == 1
                    R = obj;
                else
                    R = c * obj;
                end
            else
                R = pdmat_udiag(obj.dim, sum(obj.dv .* c));
            end
        end
        
        
        function Q = quadterm(obj, X)
            % compute the quadratic term on given samples
            
            d = obj.dim;
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == d)
                error('pdmat_udiag:invalidarg', ...
                    'X should be a numeric matrix of size d x n.');
            end
                        
            if d == 1
                Q = obj.dv' * (X.^2);
            else
                Q = obj.dv' * sum(X.^2, 1);
            end
        end

        
        function R = fullform(obj)
            % convert the matrices to full d x d x m form
            
            d = obj.dim;
            if d == 1
                R = makediag(obj.dv);
            else
                R = makediag(obj.dv(ones(d, 1), :));
            end
        end

        
        function R = chol(obj)
            % compute the (lower) Cholesky transform of all matrices
            
            R = cholmat_udiag(obj.dim, sqrt(obj.dv));
        end
        
        
        function R = sqrtm(obj)
            % compute the square root of all matrices
        
            R = pdmat_udiag(obj.dim, sqrt(obj.dv));        
        end
        
                       
        function v = logdet(obj)
            % compute the logarithm of determinants of all matrices
            
            v = obj.dim * log(obj.dv);          
        end
        
        function v = trace(obj)
            % compute the trace of all matrices
            
            v = obj.dim * obj.dv;
        end
    end
    
    
    methods(Static)
        
        function R = random(d, m)
            % Randomly generate an object
            %
            %   R = pdmat_udiag.random(d, m);
            %       randomly generates an object with m matrices of
            %       size d x d.
            %
                                                            
            R = pdmat_udiag(d, rand(1, m));
        end   
        
    end
           
end
