classdef logireg
    % The class to represent a (two-class) logistic regression problem
    % 
    %   The problem is formalized as to minimize the following 
    %   objective function:
    %
    %   sum_i w_i log (1 + exp(-y_i * theta' * x_i)) + (r2/2) * ||theta||_2^2
    %
    %   Here, both theta and x_i are d-dimensional vectors, w_i is the
    %   weight of the sample x_i, y_i is the label of x_i, whose value
    %   can be either 1 or -1.
    %
    %   Given theta, the probability that a sample x is in class with 
    %   label 1 is given by 1 / (1 + exp(-theta'*x)).
    %
    
    % Created by Dahua Lin, on Aug 7, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        dim;        % the vector space dimension (d)
        num;        % the number of samples (n)
                
        samples;    % the sample matrix [d x n]
        labels;     % the label vector [1 x n]
        weights;    % the weights of samples [1 x n]
        r2coef;     % the L2-regularization coefficient (r2)
    end
        
    methods
        
        function obj = logireg(X, y, w, r2)
            % Constructs a (two-class) logistic regression problem
            %
            %   obj = logireg(X, y, w, r2);
            %       constructs a two-class logistic regression problem
            %       with given samples and labels.
            %       
            %       Input:
            %       - X: the sample matrix of size d x n, each column is a
            %            sample
            %       - y: the label vector of size 1 x n, y(i) is the label
            %            for X(:,i).
            %       - w: the weight vector of size 1 x n. It can be input
            %            as empty, which means assigning weight 1 to 
            %            every sample.
            %       - r2:  the L2-regularization coefficient.
            %
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
                error('logireg:invalidarg', 'X should be a real matrix.');
            end
            
            [d, n] = size(X);
            
            if ~(isfloat(y) && ndims(y) == 2 && isreal(y) && ...
                    size(y, 1) == 1 && size(y, 2) == n)
                error('logireg:invalidarg', 'y should be a 1 x n real row vector.');
            end
            
            if ~(isscalar(r2) && isfloat(r2) && r2 > 0)
                error('logireg:invalidarg', 'r2 should be a positive real scalar.');
            end
            
            if ~isempty(w)
                if ~(isfloat(w) && ndims(w) == 2 && isreal(w) && ...
                        size(w, 1) == 1 && size(w, 2) == n)
                    error('logireg:invalidarg', 'w should be a 1 x n real row vector.');
                end
            end
            
            obj.dim = d;
            obj.num = n;
            obj.samples = X;
            obj.labels = y;
            obj.weights = w;
            obj.r2coef = r2;
        end
        
        
        function [v, g, H] = evaluate(obj, theta)
            % Evaluate the objective (and gradient) at a given parameter
            %
            %   v = obj.evaluate(theta);
            %       Evaluates the objective at theta.
            %
            %   [v, g] = obj.evaluate(theta);
            %       Evaluates the objective and gradient at theta.
            %
            %   [v, g, H] = obj.evaluate(obj, theta);
            %       Evaluates the objective, gradient, and Hessian matrix
            %       at theta.
            %
            
            X = obj.samples;
            y = obj.labels;
            w = obj.weights;
            r2 = obj.r2coef;
                        
            u = y .* (theta' * X);
            p = 1 ./ (1 + exp(-u));
            
            if isempty(w)
                v1 = - sum(log(p));
            else
                v1 = - log(p) * w';
            end                        
            
            v2 = (r2/2) * (theta' * theta);                        
            v = v1 + v2;
            
            if nargout >= 2
                c1 = (1 - p) .* y;
                if ~isempty(w)
                    c1 = c1 .* w;
                end
                g1 = - (X * c1');
                g2 = r2 * theta;
                
                g = g1 + g2;
            end    
            
            if nargout >= 3
                c2 = p .* (1 - p) .* (y.^2);
                if ~isempty(w)
                    c2 = c2 .* w;
                end
                H = bsxfun(@times, X, c2) * X';
                H = (H + H') / 2;
                
                d = obj.dim;
                diag_inds = 1 + (0:d-1) * (d+1);
                H(diag_inds) = H(diag_inds) + r2;
            end
        end
        
        
        function [theta, info] = solve(obj, theta0, options)
            % Solve the optimal parameter theta
            %
            %   [theta, info] = obj.solve(theta0, options);
            %       solves the optimal parameter theta using Newton's
            %       method.
            %
            %       In the input, theta0 is the initial guess of theta.
            %       options is the struct giving the optimization options.
            %       In the output, theta is the solution, and info is
            %       the procedural information of the solving process.
            %
            
            [theta, info] = newtonfmin(@(t) evaluate(obj, t), theta0, options);
        end
            
    end
    
    
    
end