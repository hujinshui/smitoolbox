classdef symat2
    % The class for representing 2x2 symmetric matrices
    %
    
    % History
    % -------
    %   - Created by Dahua Lin, on June 18, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')       
        d = 2;  % the dimension (the matrix is of size d x d)
        n;      % the number of matrices contained in the object
        
        v;      % the matrix entry values (3 x n)
    end
    
    methods
        
        %% constructor
        
        function A = symat2(M)
            % Construct an object of 2 x 2 symmetric matrices
            %
            %   A = dmat(M);
            %       constructs an object with 2 x 2 symmetric matrices.
            %
            %       For n matrix, M can be in either of the following 
            %       form:
            %       - 2 x 2 x n array with each page representing a matrix
            %       - 4 x n array with each column giving a vectorized
            %         matrix
            %       - 3 x n array with each column giving a compact
            %         vector representation [A(1,1); A(1,2); A(2,2)].
            %
            
            if isfloat(M) 
                if ndims(M) == 2
                    if size(M, 1) == 2 && size(M, 2) == 2
                        A.n = 1;
                        A.v = M([1;2;4]);
                    elseif size(M, 1) == 3
                        A.n = size(M, 2);
                        A.v = M;
                    elseif size(M, 1) == 4
                        A.n = size(M, 2);
                        A.v = M([1 2 4], :);
                    else
                        error('symat2:invalidarg', ...
                            'The size of M is invalid.');
                    end
                elseif ndims(M) == 3 && size(M, 1) == 2 && size(M, 2) == 2
                    A.n = size(M, 3);
                    V = reshape(M, 4, A.n);
                    A.v = V([1 2 4], :);
                else
                    error('symat2:invalidarg', ...
                        'The size of M is invalid.');
                end                                                
            else
                error('symat2:invalidarg', ...
                    'dv should be a numeric array.');
            end
        end
        
        
        %% Matrix retrieval
       
        function C = take(A, i)
            % Take a subset of matrix as object
            %
            %   C = A.take(i);
            %
            
            C = symat2(A.v(:, i));
        end
        
        
        function C = getm(A, i)
            % Get a particular matrix in native matrix form
            %
            %   C = A.getm(i);
            %
            
            cv = A.v(:, i);
            C = [cv(1) cv(2); cv(2) cv(3)];
        end
        
        
        function M = fullform(A)
            % Get the full form of the array of all matrices
            %
            %   M = fullform(A);
            %
            
            n_ = A.n;            
            
            if n_ == 1
                V = A.v;
                M = [V(1) V(2); V(2) V(3)];
            else
                V = A.v([1 2 2 3], :);
                M = reshape(V, [2 2 n_]);
            end
        end
        
        
        %% matrix calculation
        
        function C = plus(A, B)
            % compute the sum of two matrices
            %
            %   C = A + B;
            %
            
            if A.d == B.d
                if A.n == B.n
                    C = symat2(A.v + B.v); 
                else
                    C = symat2(bsxfun(@plus, A.v, B.v));
                end
            else
                error('MATLAB:dimagree', ...
                    'Matrix dimensions must agree.');
            end                           
        end
        
        
        function C = minus(A, B)
            % Subtract B from A
            %
            %   C = A - B;
            %
            
            if A.d == B.d
                if A.n == B.n
                    C = symat2(A.v - B.v);
                else
                    C = symat2(bsxfun(@minus, A.v, B.v));
                end
            else
                error('MATLAB:dimagree', ...
                    'Matrix dimensions must agree.');
            end
        end
                
        
        function C = times(A, B)
            % Sclalar multiplication
            %
            %   C = k .* A;
            %   C = A .* k;
            %
            %   The output is an object C with C_i = k(i) * A_i.
            %   k can be a scalar, or a vector with n elements
            %
            
            if isnumeric(A)
                k = A;
                X = B;
            else
                k = B;
                X = A;
            end
            
            v_ = X.v;
            if isscalar(k)
                C = symat2(k .* v_);
            else
                C = symat2(bsxfun(@times, k, v_));
            end
        end                
        
        
        function C = mtimes(A, B)
            % Compute matrix multiplication
            %
            %   C = A * B;
            %
            %   If both A and B are objects, then they should have
            %   the same dimension d and same number of matrices n.
            %   In this case, the output C is also an object with
            %   the same d and n, with C_i = A_i * B_i.
            %
            %   If either of A or B is a numeric matrix, then the
            %   other one should be an object with n == 1. In this
            %   case, matrix multiplication in normal sense is 
            %   performed, and the output is also a numeric matrix.
            %                                    
            
            if isnumeric(A)
                n_ = B.n;
                
                if n_ == 1
                    if size(A, 2) == 2
                        C = mtimes_sm2v2(B.v, A')';
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('symat2:singlemat', ...
                        'A single-matrix object is required.');
                end
            elseif isnumeric(B)
                n_ = A.n;
                
                if n_ == 1
                    if size(B, 1) == 2
                        C = mtimes_sm2v2(A.v, B);
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('symat2:singlemat', ...
                        'A single-matrix object is required.');
                end
            else
                n_ = A.n;
                if n_ == B.n
                    av = A.v;
                    bv = B.v;    
                    
                    if n_ == 1
                        cv = [ ...
                            av(1) * bv(1) + av(2) * bv(2); ...
                            av(1) * bv(2) + av(2) * bv(3); ...
                            av(2) * bv(2) + av(3) * bv(3)];
                        C = symat2(cv);                        
                    else
                        a1 = av(1,:);
                        a2 = av(2,:);
                        a3 = av(3,:);
                        
                        b1 = bv(1,:);
                        b2 = bv(2,:);
                        b3 = bv(3,:);
                        
                        cv = [a1 .* b1 + a2 .* b2; ...
                            a1 .* b2 + a2 .* b3; ...
                            a2 .* b2 + a3 .* b3];
                        C = symat2(cv);
                    end
                    
                else
                    error('symat2:numagree', ...
                        'A and B should have the same number of matrices.');
                end                               
            end                
        end
        
        
        function C = mldivide(A, B)
            % Compute left matrix division: inv(A) * B
            %
            %   C = A \ B
            %
            %   Here, A should be a single-matrix object,
            %   and B should be a numeric matrix with size(B,1) == A.d.
            %
                        
            if A.n == 1
                if size(B, 1) == 2
                    iv = inv2x2(A.v);
                    C = mtimes_sm2v2(iv, B);
                else
                    error('MATLAB:dimagree', ...
                        'Matrix dimensions must agree.');
                end
            else
                error('symat2:singlemat', ...
                    'A single-matrix object is required.');
            end
        end
        
        
        function Y = cmv(A, X)
            % Compute matrix-vector with corresponding vectors
            %
            %   Y = cmv(A, X);
            %
            %   Here, X should be a d x n matrix, then Y is also
            %   a d x n matrix, with Y(:,i) = A_i * X(:,i).
            %
            
            if size(X,2) == A.n
                if size(X, 1) == 2
                    Y = mtimes_sm2v2(A.v, X);
                else
                    error('MATLAB:innerdim', ...
                        'Inner matrix dimension must agree.');
                end
            else
                error('symat2:numagree', ...
                    'The number of matrices/vectors must agree.');
            end
        end
        
        
        function Y = cdv(A, X)
            % Compute left-division with corresponding vectors
            %
            %   Y = cdv(A, X);
            %
            %   Here, X should be a d x n matrix, then Y is also
            %   a d x n matrix, with Y(:,i) = A_i \ X(:,i).
            %
            
            if size(X,2) == A.n
                if size(X, 1) == 2
                    Y = mtimes_sm2v2(inv2x2(A.v), X);
                else
                    error('MATLAB:innerdim', ...
                        'Inner matrix dimension must agree.');
                end
            else
                error('symat2:numagree', ...
                    'The number of matrices/vectors must agree.');
            end
        end                                        
        
        
        %% Combine and join
        
        function C = combine(A, w)
            % Compute the weighted combination of contained matrices
            %
            %   C = combine(A);
            %   C = combine(A, w);
            %
            
            if nargin < 2 || isempty(w)
                C = symat2(sum(A.v, 2));
            else
                C = symat2(A.v * w');
            end
        end
        
        
        function C = join(A, varargin)
            % Join multiple objects together
            %
            %   C = join(A1, A2, ...)
            %
            %   A1, A2, ... should be objects with the same dimension.
            %   C is an object that contains all matrices in A1, A2, ...
            %
            
            if isempty(varargin)
                C = A;
            else
                K = numel(varargin);
                ns = zeros(1, K);
                for i = 1 : K                    
                    ns(i) = varargin{i}.n;
                end
                tn = A.n + sum(ns);
                
                vs = zeros(3, tn);
                vs(:, 1:A.n) = A.v;
                p = A.n;
                for i = 1 : K
                    vs(:, p+1 : p+ns(i)) = varargin{i}.v;
                    p = p + ns(i);
                end
                
                C = symat2(vs);
            end
        end
        
        
        %% Matrix Inverse
        
        function B = inv(A)
            % Compute the matrix inverse
            %
            %   B = inv(A)
            %
            
            B = symat2(inv2x2(A.v));
        end
        
        
        %% Characteristic numbers
        
        
        function v = lndet(A)
            % Compute the logarithm of determinant
            %
            %   v = lndet(A);
            %
            
            v = log(det2x2(A.v));
        end
        
        
        function v = trace(A)
            % Compute the trace
            %
            %   v = trace(A);
            %
            
            v = trace2x2(A.v);
        end
        
        
        %% Quadratic form
                                                
        function Q = quad(A, X, Y)
            % Compute quadratic terms 
            %
            %   Q = A.quad(X, Y);
            %
            %   It returns an n x nx matrix Q, where Q(k,i) is the
            %   x_i' * A_k * y_i, where A_k is the k-th matrix in A,
            %   and x_i is X(:,i), y_i is Y(:,i).
            %
            
            if ndims(X) == 2 && size(X,1) == 2 && ...
               ndims(Y) == 2 && size(Y,1) == 2
           
                x1 = X(1,:);
                x2 = X(2,:);
                y1 = Y(1,:);
                y2 = Y(2,:);
                
                Z = [x1 .* y1; x1 .* y2 + x2 .* y1; x2 .* y2];
                Q = A.v' * Z;               
            else
                error('symat2:invalidarg', ...
                    'X and Y should be matrices with two rows.');
            end
        end         
        
        
        %% Cholesky decomposition
        
        function Y = choltrans(A, X)
            % Apply chol(A) as transform to X
            %
            %   Y = A.choltrans(X);
            %     
            
            if A.n ~= 1
                error('symat2:invalidarg', 'A must be an single-matrix object.');
            end
            
            Lv = chol2x2(A.v);
            L = [Lv(1) 0; Lv(2) Lv(3)];            
            Y = L * X;
        end        
        
    end
    
    
    
    methods(Static)
        
        function A = randpdm(d, n)
            % Create an object with random positive definite matrix
            %
            %   A = dmat.random(d, n);
            %
            
            if d == 2
                a = rand(1, n);
                b = rand(1, n);
                t = rand(1, n) * (2*pi);
                C = mat2_by_polar(a, b, t);
                A = symat2(C);
            else
                error('symat2:randpdm:invalidarg', ...
                    'The d for symat2 can only be 2.');
            end
        end
        
        
    end
    
    
end

