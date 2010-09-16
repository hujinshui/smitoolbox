classdef gsymat
    % The class for representing a generic symmetric matrix
    %
    
    % History
    % -------
    %   - Created by Dahua Lin, on June 18, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')       
        d;      % the dimension (the matrix is of size d x d)
        n;      % the number of matrices contained in the object
        
        M;      % the matrices (d x d x n)
    end
    
    methods
        
        %% constructor
        
        function A = gsymat(M)
            % Construct an object of symmetric matrices
            %
            %   A = gsymat(dv);
            %       construct a gsymat object containing d x d symmetirc
            %       matrices. In the input, M should be an array of
            %       size d x d x n, with M(:,:,i) being the i-th
            %       symmetric matrix.
            %
            
            if isfloat(M) && ndims(M) <= 3 && size(M,1) == size(M,2)
                A.d = size(M, 1);
                A.n = size(M, 3);
                A.M = M;
            else
                error('gsymat:invalidarg', ...
                    'M should be a numeric array with ndims(M) <= 3.');
            end
        end
        
        
        %% Matrix retrieval
       
        function C = take(A, i)
            % Take a subset of matrix as object
            %
            %   C = A.take(i);
            %
            
            C = gsymat(A.M(:,:,i));
        end
        
        
        function C = getm(A, i)
            % Get a particular matrix in native matrix form
            %
            %   C = A.getm(i);
            %
            
            C = A.M(:,:,i);         
        end
        
        
        function M = fullform(A)
            % Get the full form of the array of all matrices
            %
            %   M = fullform(A);
            %
            
            M = A.M;
        end
        
        
        %% matrix calculation
        
        function C = plus(A, B)
            % compute the sum of two matrices
            %
            %   C = A + B;
            %
            
            if A.d == B.d
                if A.n == B.n
                    C = gsymat(A.M + B.M); 
                else
                    C = gsymat(bsxfun(@plus, A.M, B.M));
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
                    C = gsymat(A.M - B.M);
                else
                    C = gsymat(bsxfun(@minus, A.M, B.M));
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
            
            M_ = X.M;
            if isscalar(k)
                C = gsymat(k .* M_);
            else                
                k = reshape(k, [1, 1, length(k)]);
                C = gsymat(bsxfun(@times, k, X.M));
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
                d_ = B.d;
                
                if n_ == 1
                    if size(A, 2) == d_
                        C = A * B.M;
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('gsymat:singlemat', ...
                        'A single-matrix object is required.');
                end
            elseif isnumeric(B)
                n_ = A.n;
                d_ = A.d;
                
                if n_ == 1
                    if size(B, 1) == d_
                        C = A.M * B;
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('gsymat:singlemat', ...
                        'A single-matrix object is required.');
                end
            else
                n_ = A.n;
                d_ = A.d;
                if n_ == B.n
                    if d_ == B.d
                        if n_ == 1
                            C = gsymat(A.M * B.M);
                        else
                            AM = A.M;
                            BM = B.M;
                            CM = zeros(d_, d_, n_, class(AM(1) * BM(1)));
                            for i = 1 : n_
                                CM(:,:,i) = AM(:,:,i) * BM(:,:,i);
                            end
                            C = gsymat(CM);
                        end
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('gsymat:numagree', ...
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
            
            d_ = A.d;
            
            if A.n == 1
                if size(B, 1) == d_
                    if d_ == 1
                        C = B * (1 / A.M);
                    else
                        C = A.M \ B;
                    end
                else
                    error('MATLAB:dimagree', ...
                        'Matrix dimensions must agree.');
                end
            else
                error('gsymat:singlemat', ...
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
            
            n_ = A.n;
            if size(X, 2) == n_
                if n_ == 1
                    Y = A.M * X;
                else
                    AM = A.M;
                    Y = zeros(A.d, n_, class(AM(1) * X(1)));
                    for i = 1 : n_
                        Y(:,i) = AM(:,:,i) * X(:,i);
                    end
                end
            else
                error('gsymat:numagree', ...
                    'The numbers of matrices/vectors must agree.');
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
            
            n_ = A.n;
            if size(X, 2) == n_
                if n_ == 1
                    Y = A.M \ X;
                else
                    AM = A.M;
                    Y = zeros(A.d, n_, class(AM(1) * X(1)));
                    for i = 1 : n_
                        Y(:,i) = AM(:,:,i) \ X(:,i);
                    end
                end
            else
                error('gsymat:numagree', ...
                    'The numbers of matrices/vectors must agree.');
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
                C = gsymat(sum(A.M, 3));
            else
                d_ = A.d;
                n_ = A.n;
                AM = reshape(A.M, d_^2, n_);
                C = gsymat(reshape(AM * w', d_, d_));
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
                d_ = A.d;
                for i = 1 : K
                    if varargin{i}.d ~= d_
                        error('gsymat:invalidarg', ...
                            'All matrices should have the same dimension.');
                    end
                    ns(i) = varargin{i}.n;
                end
                tn = A.n + sum(ns);
                
                CM = zeros(d_, d_, tn);
                CM(:, :, 1:A.n) = A.M;
                p = A.n;
                for i = 1 : K
                    CM(:, :, p+1 : p+ns(i)) = varargin{i}.M;
                    p = p + ns(i);
                end
                
                C = gsymat(CM);
            end
        end
        
        
        %% Matrix Inverse
        
        function B = inv(A)
            % Compute the matrix inverse
            %
            %   B = inv(A)
            %
            
            n_ = A.n;
            if n_ == 1
                B = gsymat(inv(A.M));
            else
                d_ = A.d;
                AM = A.M;
                BM = zeros(d_, d_, n_, class(AM));
                for i = 1 : n_
                    BM(:,:,i) = inv(AM(:,:,i));
                end
                B = gsymat(BM);
            end
        end
        
        
        
        %% Characteristic numbers
        
        
        function v = lndet(A)
            % Compute the logarithm of determinant
            %
            %   v = lndet(A);
            %
            
            if A.d == 1
                v = log(reshape(A.M, 1, A.n));
            else
                n_ = A.n;
                AM = A.M;
                v = zeros(1, n_, class(AM));
                for i = 1 : n_
                    v(i) = lndet(AM(:,:,i));
                end
            end
        end
        
        
        function v = trace(A)
            % Compute the trace
            %
            %   v = trace(A);
            %
            
            if A.d == 1
                v = reshape(A.M, 1, A.n);
            else
                n_ = A.n;
                d_ = A.d;
                if n_ == 1
                    v = sum(A.M(1 + (d_+1) * (0:d_-1)));
                else
                    AM = reshape(A.M, d_ * d_, n_);
                    v = sum(AM(1 + (d_+1) * (0:d_-1), :), 1);
                end
            end
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
            
            n_ = A.n;
            m = size(X, 2);
            AM = A.M;
            
            Q = zeros(n_, m, class(AM(1) * X(1)));
            for i = 1 : n_
                Q(i, :) = sum(X .* (AM(:,:,i) * Y), 1);
            end
        end                                       
        
        
        %% Cholesky decomposition
        
        function Y = choltrans(A, X)
            % Apply chol(A) as transform to X
            %
            %   Y = A.choltrans(X);
            %     
            
            if A.n ~= 1
                error('gsymat:invalidarg', 'A must be an single-matrix object.');
            end
            
            Y = chol(A.M, 'lower') * X;
        end        
        
    end
    
    
    
    methods(Static)
        
        function A = randpdm(d, n)
            % Create an object with random positive definite matrix
            %
            %   A = gsymat.random(d, n);
            %
            
            CM = zeros(d, d, n);
            for i = 1 : n
                R = orth(randn(d, d));
                dv = rand(d, 1);
                Mi = R' * bsxfun(@times, dv, R);
                Mi = 0.5 * (Mi + Mi');
                CM(:,:,i) = Mi;
            end
            
            A = gsymat(CM);                        
        end
        
        
    end
    
    
end

