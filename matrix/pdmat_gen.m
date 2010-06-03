classdef pdmat_gen < pdmat
    % The object to represent complete Gaussian matrices
    %
    
    % Created by Dahua Lin, on Mar 22, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties
        mats;       % the matrices [d x d x m]
    end
    
    methods
        function obj = pdmat_gen(M)
            % construct an object from the given array M
            %
            %   obj = pdmat_gen(M);
            %       constructs an object that contains the Gaussian
            %       matrices given in M.
            %
            %       M should be a d x d x m array with M(:,:,i)
            %       being the i-th matrix.
            %
                        
            if ~(isfloat(M) && ndims(M) <= 3 && size(M,1) == size(M,2))
                error('pdmat_gen:invalidarg', ...
                    'M should be a either a d x d matrix or a d x d x m array.');
            end
                                    
            d = size(M, 1);
            obj = obj@pdmat(d, size(M, 3));
            
            obj.mats = M;
        end
        
        
        function R = take(obj, i)
            % takes a subset of the matrices contained in the object
            
            R = pdmat_gen(obj.mats(:,:,i));
        end

        
        function R = inv(obj)
            % compute the inverse of all matrices 
                        
            m = obj.num;
            d = obj.dim;
            M0 = obj.mats;
            
            if d == 1
                M = 1 ./ M0;
            elseif m == 1
                M = inv(M0);
            else
                d = obj.dim;
                M = zeros([d, d, m], class(obj.mats));
                for i = 1 : m
                    M(:,:,i) = inv(M0(:,:,i));
                end
            end
            
            R = pdmat_gen(M);
        end
        
        
        function b = is_subtype(obj, obj2) %#ok<MANU>
            % tests whether obj2 is of a sub-type of the class pdmat_udiag.
            
            b = isa(obj2, 'pdmat');
        end
        
        
        function Y = cmult(obj, X)
            % compute the product with corresponding vectors
            
            d = obj.dim;
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [d m]))
                error('pdmat_gen:cmult:invalidarg', ...
                    'X should be a d x m matrix.');
            end
            
            M = obj.mats;
            
            if m == 1
                Y = M * X;
            elseif d == 1
                Y = reshape(M, 1, m) .* X;
            else
                if isa(X, 'single') || isa(M, 'single')
                    Yclass = 'single';
                else
                    Yclass = 'double';
                end
                
                Y = zeros(d, m, Yclass);
                for i = 1 : m
                    Y(:,i) = M(:,:,i) * X(:,i);
                end
            end
        end        
        
        
        function Y = cldiv(obj, X)
            % compute the left matrix division with corresponding vectors
            
            d = obj.dim;
            m = obj.num;
            
            if ~(isfloat(X) && isequal(size(X), [d m]))
                error('pdmat_gen:cldiv:invalidarg', ...
                    'X should be a d x m matrix.');
            end
            
            M = obj.mats;
            
            if m == 1
                Y = M \ X;
            elseif d == 1
                Y = X ./ reshape(M, 1, m);
            else
                if isa(X, 'single') || isa(M, 'single')
                    Yclass = 'single';
                else
                    Yclass = 'double';
                end
                
                Y = zeros(d, m, Yclass);
                for i = 1 : m
                    Y(:,i) = M(:,:,i) \ X(:,i);
                end
            end            
        end
        
        
        function y = mtimes(x1, x2)
            % compute the multiplication
            
            if isa(x1, 'pdmat_gen')
                if x1.num ~= 1
                    error('pdmat_gen:mtimes:nomultmat', ...
                        'Multiplication with a multi-matrix object is not allowed.');
                end
                
                if ~(isfloat(x2) && ndims(x2) == 2 && size(x2,1) == x1.dim)
                    error('pdmat_gen:mtimes:invalidarg', ...
                        'The right multiplicand should be a d x n matrix.');
                end
                
                y = x1.mats * x2;
            else
                if ~isfloat(x1) || ...
                        ~(isscalar(x1) || x2.num == 1 || isequal(size(x1), [1 x2.num]))
                    error('pdmat_gen:mtimes:invalidarg', ...
                        'The left multiplicand should be either a scalar or a 1 x m vector.');
                end
                
                if isscalar(x1)
                    y = pdmat_gen(x1 * x2.mats);
                else
                    m = size(x1, 2);
                    M = bsxfun(@times, reshape(x1, [1 1 m]), x2.mats);
                    y = pdmat_gen(M);
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
                d = obj.dim;
                m1 = obj.num;
                m2 = rhs.num;
                
                if ~(m1 == 1 || m2 == 1 || m1 == m2)
                    error('pdmat_gen:plus:invalidarg', ...
                        'Invalid number of matrices in the objects.');
                end
                
                M0 = obj.mats;
                
                rhsclass = class(rhs);
                is_iso = strcmp(rhsclass, 'pdmat_udiag');
                is_diag = strcmp(rhsclass, 'pdmat_diag');
                
                if is_iso || is_diag
                    if d == 1
                        dv = rhs.dv;
                                                
                        M = reshape(M0, 1, m1) + dv;
                        if ~isscalar(M)
                            M = reshape(M, [1 1 numel(M)]);
                        end                            
                    else
                        if is_iso
                            dv = rhs.dv(ones(d, 1), :);
                        else
                            dv = rhs.dv;
                        end
                        
                        dins = 1 + (0:d-1) * (d+1);
                        M = reshape(M0, d*d, m1);
                        
                        if m1 == m2
                            M(dins, :) = M(dins, :) + dv;
                            M = reshape(M, [d, d, m1]);
                        elseif m1 == 1
                            M = M(:, ones(1, m2));
                            M(dins, :) = M(dins, :) + dv;
                            M = reshape(M, [d, d, m2]);
                        else % m2 == 1
                            M(dins, :) = bsxfun(@plus, M(dins, :), dv);
                            M = reshape(M, [d, d, m1]);
                        end
                    end
                else
                    if m1 == m2
                        M = M0 + fullform(rhs);
                    else
                        M = bsxfun(@plus, M0, fullform(rhs));
                    end                    
                end
                
                R = pdmat_gen(M);
            else
                error('pdmat_gen:incompatible', ...
                    'The addend is incompatible with the object.');
            end       
        end
        
        
        function R = combine(obj, c)
            % Make a single-matrix object by combining all matrices
            
            m = obj.num;    
            d = obj.dim;
            if ~(isfloat(c) && isequal(size(c), [1 m]))
                error('pdmat_gen:combine:invalidarg', ...
                    'c should be a 1 x m row vector.');
            end
                        
            if m == 1
                if c == 1
                    R = obj;
                else
                    R = c * obj;
                end
            elseif d == 1
                R = pdmat_gen(c * obj.mats(:));
            else
                d = obj.dim;
                M = reshape(obj.mats, d*d, m);
                M = M * c.'; 
                M = reshape(M, d, d);
                M = 0.5 * (M + M');
                R = pdmat_gen(M);
            end
        end
        
        
        function Y = mldivide(obj, X)
            % Compute inv(M) * X
            
            if obj.num ~= 1
                error('pdmat_gen:mldivide:momultmat', ...
                    'Do mldivide with a multi-matrix object is not allowed.');
            end
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == obj.dim)
                error('pdmat_gen:mldivide:invalidarg', ...
                    'The right multiplicand should be a d x n numeric matrix.');
            end
            
            Y = obj.mats \ X;          
        end                        
        
        
        function Q = quadterm(obj, X)
            % compute the quadratic term on given samples
            
            d = obj.dim;
            m = obj.num;
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == d)
                error('pdmat_gen:invalidarg', ...
                    'X should be a numeric matrix of size d x n.');
            end
            
            M = obj.mats;
            
            if d == 1
                X2 = X.^2;
                if m == 1
                    Q = X2 * M;
                else
                    Q = M(:) * X2;
                end
            elseif m == 1
                Q = sum(X .* (M * X), 1);
            else
                if isa(X, 'single') || isa(M, 'single')
                    Qclass = 'single';
                else
                    Qclass = 'double';
                end
                
                Q = zeros([m, size(X,2)], Qclass);
                for i = 1 : m
                    Q(i, :) = sum(X .* (M(:,:,i) * X), 1);
                end
            end            
        end

        
        function R = fullform(obj)
            % convert the matrices to full d x d x m form
            
            R = obj.mats;            
        end

        
        function R = chol(obj)
            % compute the (lower) Cholesky transform of all matrices
            
            d = obj.dim;
            m = obj.num;
            
            M = obj.mats;
            if d == 1
                T = sqrt(M);
            elseif m == 1
                T = chol(M, 'lower');
            else
                T = zeros([d, d, m], class(M));
                for i = 1 : m
                    T(:,:,i) = chol(M(:,:,i), 'lower');
                end
            end
            
            R = cholmat_gen(T);
        end
           
        
        function R = sqrtm(obj)
            % compute the square root of all matrices
            
            d = obj.dim;
            m = obj.num;
            
            M = obj.mats;
            if m == 1
                T = sqrtm(M);
            else
                T = zeros([d, d, m], class(M));
                for i = 1 : m
                    T(:,:,i) = sqrtm(M(:,:,i));
                end
            end
            
            R = pdmat_gen(T);
        end                
        
        
        function v = logdet(obj)
            % compute the logarithm of determinants of all matrices
            
            d = obj.dim;
            m = obj.num;
            M = obj.mats;
            
            if d == 1
                v = reshape(log(M), 1, m);
            elseif m == 1
                L = chol(M, 'lower');
                v = 2 * sum(log(diag(L)));
            else
                v = zeros(1, m, class(M));
                for i = 1 : m
                    L = chol(M(:,:,i), 'lower');
                    v(i) = 2 * sum(log(diag(L)));
                end
            end            
        end
        
        
        function v = trace(obj)
            % compute the trace of all matrices
            
            m = obj.num;
            d = obj.dim;
            
            M = obj.mats;
            if d == 1
                v = reshape(M, 1, m);
            elseif m == 1
                v = trace(M);
            else
                inds = 1 + (d+1) * (0:d-1);
                M = reshape(M, d*d, m);
                v = sum(M(inds, :), 1);
            end
                
        end
    end
    
    
    methods(Static)
        
        function R = random(d, m)
            % Randomly generate an object
            %
            %   R = pdmat_gen.random(d, m);
            %       randomly generates an object with m matrices of
            %       size d x d.
            %
                                                            
            M = zeros(d, d, m);
            for i = 1 : m
                es = rand(d, 1);
                U = orth(rand(d, d));
                A = U * bsxfun(@times, es, U');
                A = 0.5 * (A + A');
                M(:,:,i) = A;
            end
                
            R = pdmat_gen(M);
        end
               
    end
           
end

