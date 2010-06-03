classdef pdmat_diag < pdmat
    % The class to represent diagonal matrices        
    %
    
    % Created by Dahua Lin, on Mar 23, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties
        dv;       % the value of diagonal entries [d x m]
    end
    
    methods
        function obj = pdmat_diag(dv)
            % construct an object from the given diagonal values
            %
            %   obj = pdmat_diag(dv);
            %       constructs an object that contains the d x d diagonal
            %       matrices with diagonal elements given by dv.
            %
            %       dv should be a d x m row vector, then the i-th            
            %       matrix is diag(dv(:,i)).
            %
                                 
            if ~(isfloat(dv) && ndims(dv) == 2)
                error('pdmat_diag:invalidarg', 'dv should be numeric matrix.');
            end
                               
            d = size(dv, 1);
            obj = obj@pdmat(d, size(dv, 2));            
            obj.dv = dv;
        end
        
        
        function R = take(obj, i)
            % takes a subset of the matrices contained in the object
            
            R = pdmat_diag(obj.dv(:, i));
        end

        
        function R = inv(obj)
            % compute the inverse of all matrices 
                                                
            R = pdmat_diag(1 ./ obj.dv);
        end
        
        
        function b = is_subtype(obj, obj2) %#ok<MANU>
            % tests whether obj2 is of a sub-type of the class pdmat_udiag.
            
            b = isa(obj2, 'pdmat_udiag') || isa(obj2, 'pdmat_diag');
        end
        
        
        function Y = cmult(obj, X)
            % compute the product with corresponding vectors
            
            d = obj.dim;
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [d m]))
                error('pdmat_diag:cmult:invalidarg', ...
                    'X should be a d x m matrix.');
            end
                                    
            Y = obj.dv .* X;
        end        
        
        
        function Y = cldiv(obj, X)
            % compute the left matrix division with corresponding vectors
            
            d = obj.dim;
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [d m]))
                error('pdmat_diag:cldiv:invalidarg', ...
                    'X should be a d x m matrix.');
            end
                                    
            Y = X ./ obj.dv;
        end
        
        
        function y = mtimes(x1, x2)
            % compute the multiplication
            
            if isa(x1, 'pdmat_diag')            
                if x1.num ~= 1
                    error('pdmat_diag:mtimes:nomultmat', ...
                        'Multiplication with a multi-matrix object is not allowed.');
                end
            
                if ~(isfloat(x2) && ndims(x2) == 2 && size(x2,1) == x1.dim)
                    error('pdmat_diag:mtimes:invalidarg', ...
                        'The right multiplicand should be a d x n matrix.');
                end
                
                if x1.dim == 1
                    y = x1.dv * x2;
                else
                    y = bsxfun(@times, x1.dv, x2);
                end
            else
                if ~isfloat(x1) || ...
                        ~(isscalar(x1) || x2.num == 1 || isequal(size(x1), [1 x2.num]))
                    error('pdmat_diag:mtimes:invalidarg', ...
                        'The left multiplicand should be either a scalar or a 1 x m vector.');
                end
                
                if isscalar(x1)
                    y = pdmat_diag(x1 * x2.dv);
                else
                    y = pdmat_diag(bsxfun(@times, x1, x2.dv));
                end
            end
        end       
        
        
        function R = plus(obj, rhs)
            % compute the addition of two matrices
                
            if obj.dim ~= rhs.dim
                error('pdmat_gen:plus:invalidarg', ...
                    'The dimensions of addends are not the same.');
            end
            
            if is_subtype(obj, rhs)
                if ~(obj.num == 1 || rhs.num == 1 || obj.num == rhs.num)
                    error('pdmat_diag:plus:invalidarg', ...
                        'Invalid number of matrices in the objects.');
                end
                
                odv = obj.dv;
                rdv = rhs.dv;
                
                if isscalar(odv) || isscalar(rdv) || all(size(odv) == size(rdv))                        
                    adv = odv + rdv;
                else
                    adv = bsxfun(@plus, odv, rdv);
                end
                    
                R = pdmat_diag(adv);
                
            elseif is_subtype(rhs, obj)                
                R = plus(rhs, obj);
                
            else
                error('pdmat_diag:incompatible', ...
                    'The addend is incompatible with the object.');
            end           
        end
        
        
        function Y = mldivide(obj, X)
            % Compute inv(M) * X
            
            if obj.num ~= 1
                error('pdmat_diag:mldivide:momultmat', ...
                    'Do mldivide with a multi-matrix object is not allowed.');
            end
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == obj.dim)
                error('pdmat_diag:mldivide:invalidarg', ...
                    'The right multiplicand should be a d x n numeric matrix.');
            end
            
            Y = bsxfun(@times, X, 1 ./ obj.dv);          
        end
        
                                
        function R = combine(obj, c)
            % Make a single-matrix object by combining all matrices
            
            m = obj.num;            
            if ~(isfloat(c) && isequal(size(c), [1 m]))
                error('pdmat_diag:combine:invalidarg', ...
                    'c should be a 1 x m row vector.');
            end
            
            if m == 1 && c == 1
                R = obj;
            else
                R = pdmat_diag(obj.dv * c.');
            end
        end
        
                                    
        function Q = quadterm(obj, X)
            % compute the quadratic term on given samples
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == obj.dim)
                error('pdmat_diag:invalidarg', ...
                    'X should be a numeric matrix of size d x n.');
            end
                        
            Q = obj.dv.' * (X.^2);
        end

        
        function R = fullform(obj)
            % convert the matrices to full d x d x m form
            
            R = makediag(obj.dv);
        end

        
        function R = chol(obj)
            % compute the (lower) Cholesky transform of all matrices
            
            R = cholmat_diag(sqrt(obj.dv));
        end
                        
        
        function R = sqrtm(obj)
            % compute the square root of all matrices
        
            R = pdmat_diag(sqrt(obj.dv));        
        end
        
        
        function v = logdet(obj)
            % compute the logarithm of determinants of all matrices
        
            if obj.dim == 1
                v = log(obj.dv);
            else
                v = sum(log(obj.dv), 1);
            end
        end
        
        
        function v = trace(obj)
            % compute the trace of all matrices
            
            if obj.dim == 1
                v = obj.dv;
            else
                v = sum(obj.dv, 1);
            end            
        end        
        
    end
    
    
    methods(Static)
        
        function R = random(d, m)
            % Randomly generate an object
            %
            %   R = pdmat_diag.random(d, m);
            %       randomly generates an object with m matrices of
            %       size d x d.
            %
                                                            
            R = pdmat_diag(rand(d, m));
        end   
        
    end
           
end
