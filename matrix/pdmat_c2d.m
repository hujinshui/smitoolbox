classdef pdmat_c2d < pdmat
    % The object to represent 2x2 positive definite matrices 
    %
    % Each matrix is represented by three values [c11, c12, c22]  
    %
    % This class is specially tailored for 2x2 matrices with efficient
    % implementation
    %    
    
    % Created by Dahua Lin, on Mar 23, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties
        v;       % the values [3 x m], each column in form of [c11; c12; c22]
    end
    
    methods
        function obj = pdmat_c2d(v)
            % construct an object from the given entry values
            %
            %   obj = pdmat_c2d(v);
            %       constructs an object that contains the 2 x 2 Gaussian
            %       matrices with entry values given by v.
            %
            %       v should be a 3 x m row vector, then the i-th            
            %       matrix is [v(1,i) v(2,i); v(2,i) v(3,i)].
            %
                                 
            if ~(isfloat(v) && ndims(v) == 2 && size(v,1) == 3)
                error('pdmat_c2d:invalidarg', ...
                    'v should be 3 x m numeric matrix.');
            end
                                           
            obj = obj@pdmat(2, size(v, 2));            
            obj.v = v;
        end
        
        
        function R = take(obj, i)
            % takes a subset of the matrices contained in the object
            
            R = pdmat_c2d(obj.v(:, i));
        end

        
        function R = inv(obj)
            % compute the inverse of all matrices 
                                    
            v0 = obj.v;
            dv = v0(1,:) .* v0(3,:) - v0(2,:).^2;
            iv = bsxfun(@times, [v0(3,:); -v0(2,:); v0(1,:)], 1./ dv);
            
            R = pdmat_c2d(iv);
        end
        
        
        function b = is_subtype(obj, obj2) %#ok<MANU>
            % tests whether obj2 is of a sub-type of the class pdmat_udiag.
            
            b = isa(obj2, 'pdmat') && obj2.dim == 2;
        end
        
        
        function Y = cmult(obj, X)
            % compute the product with corresponding vectors
            
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [2 m]))
                error('pdmat_c2d:cmult:invalidarg', ...
                    'X should be a 2 x m matrix.');
            end
                                    
            V = obj.v;
            y1 = sum(V([1 2], :) .* X, 1);
            y2 = sum(V([2 3], :) .* X, 1);
            
            Y = [y1; y2];            
        end        
        
        
        function Y = cldiv(obj, X)
            % compute the left matrix division with corresponding vectors
            
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [2 m]))
                error('pdmat_c2d:cldiv:invalidarg', ...
                    'X should be a 2 x m matrix.');
            end
                                    
            V = obj.v;
            a = V(1, :);
            b = V(2, :);
            c = V(3, :);
            k = 1 ./ (a .* c - b.^2);
            
            x1 = X(1, :);
            x2 = X(2, :);
            
            y1 = k .* (c .* x1 - b .* x2);
            y2 = k .* (a .* x2 - b .* x1);
                        
            Y = [y1; y2];         
        end
        
        
        function y = mtimes(x1, x2)
            % compute the multiplication
            
            if isa(x1, 'pdmat_c2d')            
                if x1.num ~= 1
                    error('pdmat_c2d:mtimes:nomultmat', ...
                        'Multiplication with a multi-matrix object is not allowed.');
                end
            
                if ~(isfloat(x2) && ndims(x2) == 2 && size(x2,1) == 2)
                    error('pdmat_c2d:mtimes:invalidarg', ...
                        'The right multiplicand should be a 2 x n matrix.');
                end
                
                r1 = x2(1, :);
                r2 = x2(2, :);
                v0 = x1.v;
                
                y = [v0(1,:) .* r1 + v0(2,:) .* r2; v0(2,:) .* r1 + v0(3,:) .* r2];
            else
                if ~isfloat(x1) || ~(isscalar(x1) || x2.num == 1 || isequal(size(x1), [1 x2.num]))
                    error('pdmat_c2d:mtimes:invalidarg', ...
                        'The left multiplicand should be either a scalar or a 1 x m vector.');
                end
                
                if isscalar(x1)
                    y = pdmat_c2d(x1 * x2.v);
                else
                    y = pdmat_c2d(bsxfun(@times, x1, x2.v));
                end
            end
        end       
        
        
        function R = plus(obj, rhs)
            % compute the addition of two matrices
                                                
            if is_subtype(obj, rhs)
                m1 = obj.num;
                m2 = rhs.num;
                
                if ~(m1 == 1 || m2 == 1 || m1 == m2)
                    error('pdmat_c2d:plus:invalidarg',  ...
                        'Invalid number of matrices in the objects.');
                end
                
                v0 = obj.v;
                
                switch class(rhs)
                    case 'pdmat_c2d'
                        if m1 == m2
                            va = v0 + rhs.v;
                        else
                            va = bsxfun(@plus, v0, rhs.v);
                        end
                        
                    case 'pdmat_udiag'
                        va = v0;
                        dv2 = rhs.dv;
                        if m1 < m2  % m1 == 1 and m2 > 1
                            va = va(:, ones(1, m2));
                        end
                        va(1, :) = va(1, :) + dv2;
                        va(3, :) = va(3, :) + dv2;
                        
                    case 'pdmat_diag'                        
                        dv2 = rhs.dv;
                        v2 = [dv2(1,:); zeros(1, m2); dv2(2,:)];
                        if m1 == m2
                            va = v0 + v2;
                        else
                            va = bsxfun(@plus, v0, v2);
                        end
                        
                    otherwise
                        v2 = reshape(fullform(rhs), 4, m2);
                        v2 = v2([1 2 4], :);
                        if m1 == m2
                            va = v0 + v2;
                        else
                            va = bsxfun(@plus, v0, v2);
                        end
                end
                
                R = pdmat_c2d(va);
            else
                error('pdmat_c2d:incompatible', ...
                    'The addend is incompatible with the object.');
            end            
        end
        
        
        function Y = mldivide(obj, X)
            % Compute inv(M) * X
            
            if obj.num ~= 1
                error('pdmat_c2d:mldivide:momultmat', ...
                    'Do mldivide with a multi-matrix object is not allowed.');
            end
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == 2)
                error('pdmat_c2d:mldivide:invalidarg', ...
                    'The right multiplicand should be a 2 x n matrix.');
            end
            
            v0 = obj.v;
            a = v0(1);
            b = v0(2);
            c = v0(3);
            
            k = 1 / (a * c - b^2);
            
            x1 = X(1, :);
            x2 = X(2, :);
            
            Y = [k * (c * x1 - b * x2); k * (a* x2 - b * x1)]; 
        end
        
        
        
        function R = combine(obj, c)
            % Make a single-matrix object by combining all matrices
            
            m = obj.num;            
            if ~(isfloat(c) && isequal(size(c), [1 m]))
                error('pdmat_c2d:combine:invalidarg', ...
                    'c should be a 1 x m row vector.');
            end
            
            if m == 1 && c == 1
                R = obj;
            else
                R = pdmat_c2d(obj.v * c.');
            end
        end
        
                      
        function Q = quadterm(obj, X)
            % compute the quadratic term on given samples
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == 2)
                error('pdmat_c2d:invalidarg', ...
                    'X should be a numeric matrix of size 2 x n.');
            end
                        
            A = [X(1,:).^2; 2 * X(1,:) .* X(2,:); X(2,:).^2];
            Q = obj.v.' * A;
        end

        
        function R = fullform(obj)
            % convert the matrices to full d x d x m form
                        
            m = obj.num;            
            V = obj.v;
            
            if m == 1
                R = [V(1), V(2); V(2), V(3)];
            else
                d = obj.dim;
                R = reshape(V([1 2 2 3], :), [d d m]);
            end
        end

        
        function R = chol(obj)
            % compute the (lower) Cholesky transform of all matrices
            
            C = chol2x2(obj.v);
                                    
            m = obj.num;
            if m == 1
                T = [C(1) 0; C(2) C(3)];
            else
                T = [C([1 2], :); zeros(1, m); C(3,:)];
                T = reshape(T, [2 2 m]);
            end
            R = cholmat_gen(T);
        end
        
        
        function R = sqrtm(obj)
            % compute the square root of all matrices
        
            R = pdmat_c2d(sqrtm2x2(obj.v));        
        end
                
        
        function dv = logdet(obj)
            % compute the logarithm of determinants of all matrices
                        
            dv = log(det2x2(obj.v));
        end
        
        function tv = trace(obj)
            % Compute the trace of all matrices
            
            v0 = obj.v;
            tv = v0(1,:) + v0(3,:);
        end
    end
    
    
    methods(Static)
        
        function R = random(d, m)
            % Randomly generate an object
            %
            %   R = pdmat_c2d.random(d, m);
            %       randomly generates an object with m matrices of
            %       size d x d.
            %
            
            if ~isequal(d, 2)
                error('pdmat_c2d:random:invalidarg', ...
                    'd can only be 2 for pdmat_c2d.');
            end
                                                            
            h = rand(3, m);
            a = h(1, :); 
            b = h(2, :);
            rou = 2 * h(3, :) - 1;
            V = [a .^ 2; a .* b .* rou; b .^ 2];
            R = pdmat_c2d(V);
        end   
        
    end
                      
end
