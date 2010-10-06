classdef dmat
    % The class for representing a diagonal matrix
    %
    
    % History
    % -------
    %   - Created by Dahua Lin, on June 18, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')       
        d;      % the dimension (the matrix is of size d x d)
        n;      % the number of matrices contained in the object
        
        dv;     % the diagonal value (d x n)
    end
    
    methods
        
        %% constructor
        
        function A = dmat(dv)
            % Construct an object of diagonal matrices
            %
            %   A = dmat(dv);
            %       construct a dmat object containing d x d diagonal 
            %       matrices, where the diagonal values of the i-th 
            %       diagonal matrix is given by dv(i), and the number
            %       of diagonal matrices is size(dv, 2).                        
            %
            
            if isfloat(dv) && ndims(dv) == 2                                    
                [A.d, A.n] = size(dv);
                A.dv = dv;
            else
                error('dmat:invalidarg', ...
                    'dv should be a numeric matrix.');
            end
        end
        
        
        %% Matrix retrieval
       
        function C = take(A, i)
            % Take a subset of matrix as object
            %
            %   C = A.take(i);
            %
            
            C = dmat(A.dv(:, i));
        end
        
        
        function C = getm(A, i)
            % Get a particular matrix in native matrix form
            %
            %   C = A.getm(i);
            %
            
            d_ = A.d;
            dv_ = A.dv;
            C = zeros(d_, class(dv_));
            C(1 + (0:d_-1) * (d_+1)) = dv_(:,i);            
        end
        
        
        function M = fullform(A)
            % Get the full form of the array of all matrices
            %
            %   M = fullform(A);
            %
            
            d_ = A.d;                        
            n_ = A.n;            
            
            if n_ == 1                            
                M = zeros(d_, d_, class(A.dv));
                M(1 + (0:d_-1) * (d_+1)) = A.dv;
            elseif d_ == 1
                M = reshape(A.dv, [1, 1, n_]);
            else
                M = zeros(d_ * d_, n_, class(A.dv));
                M(1 + (0:d_ - 1) * (d_ + 1), :) = A.dv;
                M = reshape(M, [d_, d_, n_]);
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
                    C = dmat(A.dv + B.dv); 
                else
                    C = dmat(bsxfun(@plus, A.dv, B.dv));
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
                    C = dmat(A.dv - B.dv);
                else
                    C = dmat(bsxfun(@minus, A.dv, B.dv));
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
            
            dv_ = X.dv;
            if isscalar(k) || X.d == 1
                C = dmat(k .* dv_);
            else
                C = dmat(bsxfun(@times, k, dv_));
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
                        if d_ == 1
                            C = A * B.dv;                            
                        else
                            C = bsxfun(@times, A, B.dv.');
                        end
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('dmat:singlemat', ...
                        'A single-matrix object is required.');
                end
            elseif isnumeric(B)
                n_ = A.n;
                d_ = A.d;
                
                if n_ == 1
                    if size(B, 1) == d_
                        if d_ == 1
                            C = A.dv * B;
                        else
                            C = bsxfun(@times, A.dv, B);
                        end
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('dmat:singlemat', ...
                        'A single-matrix object is required.');
                end
            else
                d_ = A.d;
                if A.n == B.n
                    if d_ == B.d
                        C = dmat(A.dv .* B.dv);
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('dmat:numagree', ...
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
                        C = B * (1 / A.dv);
                    else
                        C = bsxfun(@times, B, 1 ./ A.dv);
                    end
                else
                    error('MATLAB:dimagree', ...
                        'Matrix dimensions must agree.');
                end
            else
                error('dmat:singlemat', ...
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
            
            Y = X .* A.dv;
        end
        
        
        function Y = cdv(A, X)
            % Compute left-division with corresponding vectors
            %
            %   Y = cdv(A, X);
            %
            %   Here, X should be a d x n matrix, then Y is also
            %   a d x n matrix, with Y(:,i) = A_i \ X(:,i).
            %
            
            Y = X ./ A.dv;
        end                                        
        
        
        %% Combine and join
        
        function C = combine(A, w)
            % Compute the weighted combination of contained matrices
            %
            %   C = combine(A);
            %   C = combine(A, w);
            %
            
            if nargin < 2 || isempty(w)
                C = dmat(sum(A.dv, 2));
            else
                C = dmat(A.dv * w');
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
                        error('dmat:invalidarg', ...
                            'All matrices should have the same dimension.');
                    end
                    ns(i) = varargin{i}.n;
                end
                tn = A.n + sum(ns);
                
                dvs = zeros(d_, tn);
                dvs(:, 1:A.n) = A.dv;
                p = A.n;
                for i = 1 : K
                    dvs(:, p+1 : p+ns(i)) = varargin{i}.dv;
                    p = p + ns(i);
                end
                
                C = dmat(dvs);
            end
        end
        
        
        %% Matrix Inverse
        
        function B = inv(A)
            % Compute the matrix inverse
            %
            %   B = inv(A)
            %
            
            B = dmat(1 ./ A.dv);
        end
        
        
        
        %% Characteristic numbers
        
        
        function v = lndet(A)
            % Compute the logarithm of determinant
            %
            %   v = lndet(A);
            %
            
            if A.d == 1
                v = log(A.dv);
            else
                v = sum(log(A.dv), 1);
            end
        end
        
        
        function v = trace(A)
            % Compute the trace
            %
            %   v = trace(A);
            %
            
            if A.d == 1
                v = A.dv;
            else
                v = sum(A.dv, 1);
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
            
            Q = A.dv' * (X .* Y);       
        end    
        
        
        %% Cholesky decomposition
        
        function Y = choltrans(A, X)
            % Apply chol(A) as transform to X
            %
            %   Y = A.choltrans(X);
            %     
            
            if A.n ~= 1
                error('dmat:invalidarg', 'A must be an single-matrix object.');
            end
            
            Y = bsxfun(@times, sqrt(A.dv), X);
        end        
        
    end
    
    
    
    methods(Static)
        
        function A = randpdm(d, n)
            % Create an object with random positive definite matrix
            %
            %   A = dmat.random(d, n);
            %
            
            A = dmat(rand(d, n));                        
        end
        
        
    end
    
    
end

