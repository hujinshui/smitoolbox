classdef udmat
    % The class for representing a uniform diagonal matrix (i.e. dv * I)
    %
    
    % History
    % -------
    %   - Created by Dahua Lin, on June 11, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')       
        d;      % the dimension (the matrix is of size d x d)
        n;      % the number of matrices contained in the object
        
        dv;     % the diagonal value (1 x n)
    end
    
    methods
        
        %% constructor
        
        function A = udmat(d, dv)
            % Construct an object of uniform diagonal matrices
            %
            %   A = udmat(d, dv);
            %       construct a udmat object containing d x d
            %       uniformly diagonal matrices, where the
            %       diagonal value of the i-th diagonal matrix
            %       is given by dv(i). 
            %
            %       dv should be an 1 x n vector.
            %
            
            if isfloat(dv) && ndims(dv) == 2 && size(dv,1) == 1            
                A.d = d;
                A.n = size(dv, 2);
                A.dv = dv;
            else
                error('udmat:invalidarg', 'dv should be a row vector.');
            end
        end
        
        
        %% Matrix retrieval
       
        function C = take(A, i)
            % Take a subset of matrix as object
            %
            %   C = A.take(i);
            %
            
            C = udmat(A.d, A.dv(i));
        end
        
        
        function C = getm(A, i)
            % Get a particular matrix in native matrix form
            %
            %   C = A.getm(i);
            %
            
            dim = A.d;
            C = zeros(dim, class(A.dv));
            C(1 + (0:dim-1) * (dim+1)) = A.dv(i);            
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
                M(1 + (0:d_ - 1) * (d_ + 1), :) = A.dv(ones(d_, 1), :);
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
                C = udmat(A.d, A.dv + B.dv);                 
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
                C = udmat(A.d, A.dv - B.dv);                 
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
            
            C = udmat(X.d, k .* X.dv);                            
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
                if B.n == 1
                    if size(A, 2) == B.d
                        C = B.dv * A;
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('udmat:singlemat', ...
                        'A single-matrix object is required.');
                end
            elseif isnumeric(B)
                if A.n == 1
                    if size(B, 1) == A.d
                        C = A.dv * B;
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('udmat:singlemat', ...
                        'A single-matrix object is required.');
                end
            else
                d_ = A.d;
                if A.n == B.n
                    if d_ == B.d
                        C = udmat(d_, A.dv .* B.dv);
                    else
                        error('MATLAB:innerdim', ...
                            'Inner matrix dimension must agree.');
                    end
                else
                    error('udmat:numagree', ...
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
                if A.d == size(B, 1)
                    C = B * (1 / A.dv);
                else
                    error('MATLAB:dimagree', ...
                        'Matrix dimensions must agree.');
                end
            else
                error('udmat:singlemat', ...
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
            if n_ == size(X, 2)
                if A.d == size(X, 1)
                    Y = bsxfun(@times, A.dv, X);
                else
                    error('MATLAB:dimagree', ...
                        'Matrix dimensions must agree.');
                end
            else
                error('udmat:numagree', ...
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
            if n_ == size(X, 2)
                if A.d == size(X, 1)
                    Y = bsxfun(@times, 1 ./ A.dv, X);
                else
                    error('MATLAB:dimagree', ...
                        'Matrix dimensions must agree.');
                end
            else
                error('udmat:numagree', ...
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
                C = udmat(A.d, sum(A.dv));
            else
                C = udmat(A.d, A.dv * w');
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
                        error('udmat:invalidarg', ...
                            'All matrices should have the same dimension.');
                    end
                    ns(i) = varargin{i}.n;
                end
                tn = A.n + sum(ns);
                
                dvs = zeros(1, tn);
                dvs(1:A.n) = A.dv;
                p = A.n;
                for i = 1 : K
                    dvs(p+1 : p+ns(i)) = varargin{i}.dv;
                    p = p + ns(i);
                end
                
                C = udmat(d_, dvs);
            end
        end
        
        
        %% Matrix Inverse
        
        function B = inv(A)
            % Compute the matrix inverse
            %
            %   B = inv(A)
            %
            
            B = udmat(A.d, 1 ./ A.dv);
        end
        
        
        
        
        %% Characteristic numbers
        
        
        function v = lndet(A)
            % Compute the logarithm of determinant
            %
            %   v = lndet(A);
            %
            
            v = log(A.dv) * A.d;
        end
        
        
        function v = trace(A)
            % Compute the trace
            %
            %   v = trace(A);
            %
            
            v = A.dv * A.d;
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
            Q = sum(X .* Y, 1);
            if n_ == 1
                Q = Q * A.dv;
            else
                Q = bsxfun(@times, Q, A.dv.');
            end            
        end                                       
        
    end
    
    
    
    methods(Static)
        
        function A = randpdm(d, n)
            % Create an object with random positive definite matrix
            %
            %   A = udmat.random(d, n);
            %
            
            A = udmat(d, rand(1, n));                        
        end
        
        
    end
    
    
end

