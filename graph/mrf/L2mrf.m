classdef L2mrf
    % The class to represent Markov Random Field based on L2-potentials
    %
    %   The energy of an L2-MRF with affinity matrix W is 
    %   defined to be
    %
    %       E(X) = (1/2) * sum_{i<j} W(i,j) ||x_i - x_j||_2^2
    %            = (1/2) * trace(X * L * X');
    %
    %   Here, X = [x_1, ..., x_n], and L is the Laplacian matrix.
    %
    %   The energy E can be converted to pdf value as exp(-E)/Z. 
    %
    
    % Created by Dahua Lin, on Sep 20, 2010
    %
    
    
    properties(GetAccess='public', SetAccess='private')
        nnodes; % the number of nodes (n)
        
        sw;     % the sum of weights for each node [1 x n]
        W;      % the affinity matrix (symmetric) [n x n]
        L;      % the Laplacian matrix [n x n]
    end
        
    methods
        
        function obj = L2mrf(W)
            % Construct a Gaussian MRF
            %
            %   obj = L2mrf(W);
            %       constructs an L2 MRF with its affinity matrix W.
            %       
            %       Note that W should be a symmetric matrix. W can be
            %       either a full matrix or a sparse matrix.
            %
            
            if ~(isfloat(W) && ndims(W) == 2 && isreal(W) && size(W,1) == size(W,2))
                error('gmrf:invalidarg', 'W should be a symmetric real matrix.');
            end
            
            obj.nnodes = size(W,1);
            obj.sw = full(sum(W,1));
            obj.W = W;
            obj.L = laplacemat(W);                
        end
        
        
        function e = energy(obj, X)
            % Evaluate the energy on given samples
            %
            %   e = obj.energy(X);
            %
            
            XL = X * obj.L;
            e = sum(dot(X, XL)) / 2;
        end
        
        
        function X = solve(obj, H, a, fsolver)
            % Solve the values at all nodes
            %
            %   X = obj.solve(H, a);
            %       solves the problem of minimizing the following
            %       objective:
            %
            %           E(X) - <H, X> + sum_i a(i)/2 * ||x_i||^2
            %       
            %       Here, <H, X> = sum_i h_i^T * x_i. There is a
            %       closed-form solution, given by            
            %
            %           X = (L + diag(a)) \ H
            %
            %       In the input, H should be a matrix with n columns,
            %       and a can be either a scalar or a vector of length n.
            %
            %   X = obj.solve(H, a, fsolver);
            %       specifies the solver to solve the linear equation.
            %       
            
            n = obj.nnodes;            
            if nargin < 4
                fsolver = [];
            end
            
            if ~(isfloat(H) && ndims(H) == 2 && size(H,2) == n)
                error('L2mrf:solve:invalidarg', ...
                    'H should be a numeric matrix with n columns.');
            end
            
            Q = adddiag(obj.L, a);
            
            if isempty(fsolver)
                X = H / Q;
            else
                X = fsolver(Q, H')';
            end
            
        end
        
        
        function X = smooth(obj, Y, a, fsolver)
            % Solve the problem of smoothing signals
            %
            %   X = obj.smooth(obj, Y, a);
            %   X = obj.smooth(obj, Y, a, fsolver);
            %       solves the problem of minimizing the following
            %       objective:
            %
            %           E(X) + sum_i a(i) / 2 * ||x_i - y_i||^2
            %
            %       In the input, Y should be a matrix with n columns, 
            %       and a can be either a scalar, or a row vector of 
            %       length n.
            %
            
            n = obj.nnodes;            
            if nargin < 4
                fsolver = [];
            end
            
            if ~(isfloat(Y) && ndims(Y) == 2 && size(Y,2) == n)
                error('L2mrf:solve:invalidarg', ...
                    'H should be a numeric matrix with n columns.');
            end
            
            % compute H
            
            if isscalar(a)
                if a == 1
                    H = Y;
                else
                    H = a * Y;
                end
            else
                H = bsxfun(@times, Y, a);
            end
            
            X = obj.solve(H, a, fsolver);
            
        end
    end
    
end

