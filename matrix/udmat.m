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
            %   For efficiency, argument checking is not performed
            %   by the function. It is the caller's responsibility
            %   to ensure the correctness of input arguments.
            %
            
            A.d = d;
            A.dv = dv;
            A.n = numel(dv);
        end
       
        %% operator overloading
        
        function C = plus(A, B)
            % compute the sum of two matrices
            %
            %   C = A + B;
            %
            
            if (A.d ~= B.d)
                error('udmat:invalidarg', ...
                    'A and B should have the same dimension.');
            end            
            C = udmat(A.d, A.dv + B.dv);
        end
        
        
        function C = minus(A, B)
            % Subtract B from A
            %
            %   C = A - B;
            %
            
            if (A.d ~= B.d)
                error('udmat:invalidarg', ...
                    'A and B should have the same dimension.');
            end            
            C = udmat(A.d, A.dv - B.dv);
        end
        
        
        function C = mtimes(A, B)
            % Compute scalar multiplication or matrix multiplication
            %
            %   C = alpha * A;
            %   C = A * alpha;
            %   C = A * B;
            %
            
            if isnumeric(A)  % scalar mult. with scalar on the left                
                C = udmat(B.d, A .* B.dv);                                
            elseif isnumeric(B)  % scalar mult. with scalar on the right
                C = udmat(A.d, A.dv .* B);
            else  % matrix mult.
                if (A.d ~= B.d)
                    error('udmat:invalidarg', ...
                        'A and B should have the same dimension.');
                end                
                C = udmat(A.d, A.dv .* B.dv);
            end
        end
        
        
        function C = mldivide(A, B)
            % Compute left matrix division: inv(A) * B
            %
            %   C = A \ B
            %
            
            if (A.d ~= B.d)
                error('udmat:invalidarg', ...
                    'A and B should have the same dimension.');
            end
            C = udmat(A.d, B.dv ./ A.dv);            
        end
        
        
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
        
        
        
    end
    
    
end