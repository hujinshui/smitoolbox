classdef gausslin < gcondmodel
    % The class to represent Gaussian linear model 
    %
    %   A Gaussian linear model is a conditional likelihood model
    %   defined by    
    %       p(y | x; A, mu, sigma) = N(y - Ax | mu, sigma);
    %
    %   Here N denotes a Gaussian distribution. A is the transform
    %   matrix.
    %
    %   Note that the dimension of x and that of y can be different.
    %
    
    % Created by Dahua Lin, on Mar 24, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        A;          % the transform matrix 
                    % it can be either a d x r matrix, a scalar, or empty.
                    % when A is empty, it represents identity transform.
        
        xdim;       % the (input) dimension (r);
        ydim;       % the (output) dimension (d);
        
        mu;         % the mean vector of the internal Gaussian model
        icov;       % the inverse covariance of the internal Gaussian model
        
        Cx;         % the matrix A' * inv(sigma) * A
        Rx;         % the matrix inv(sigma) * A
    end
    
    methods        
        
        function obj = gausslin(A, g)
            % Constructs a Gaussian linear model
            %
            %   obj = gausslin(A, g);
            %       constructs a Gaussian linear model with the
            %       given transform matrix and then conditional 
            %       Gaussian model g.
            %
            %       Here, A can be in either of the following forms:
            %       - a d x r matrix (r == g.dim)
            %       - a scalar (in this case, the transform matrix is A *
            %         eye(d)), here d can be determined from g.
            %
            %   `   g should be a gaussd_cp object.
            %
            
            % process input Gaussian model
            
            if isa(g, 'gaussd_cp')
                if ~(g.nmodels == 1)
                    error('gausslin:invalidarg', ...
                        'g should have a single model in it.');
                end
                gt1 = g.theta1;
                gt2 = g.theta2;
                
                if all(gt1 == 0)
                    obj.mu = 0;
                else
                    obj.mu = gt2 \ gt1;
                end
                
                obj.icov = gt2;
                
            elseif isa(g, 'gaussd_mp')
                if ~(g.nmodels == 1)
                    error('gausslin:invalidarg', ...
                        'g should have a single model in it.');
                end
                
                obj.mu = g.mu;
                gt2 = inv(g.sigma);
                obj.icov = gt2;
            else
                error('gausslin:invalidarg', ...
                    'g should be a gaussd_cp or gaussd_mp object.');
            end
                        
            % process transform A
            
            d = g.dim;
            if isempty(A)
                A = [];
                r = d;
            else
                if ~isfloat(A) || ~isreal(A) || ...
                    ~(isscalar(A) || (ndims(A) == 2 && size(A,1) == d)), ...
                    error('gausslin:invalidarg', ...
                        'A is in an invalid form.');
                end
                
                if isscalar(A)
                    r = d;obj.check_xy_args(X, Y);  
                else
                    r = size(A, 2);
                end
            end
            
            obj.A = A;
            obj.xdim = r;
            obj.ydim = d;
            obj.gauss = g;              
                       
            % pre-compute Rx and Cx
            
            if isempty(A)
                Rx_ = gt2;
                Cx_ = Rx_;
            elseif isscalar(A)
                Rx_ = A * gt2;
                Cx_ = (A^2) * gt2;
            else
                Rx_ = gt2 * A;
                Cx_ = A' * Rx_;
                Cx_ = 0.5 * (Cx_ + Cx_');
            end
                
            obj.Rx_ = Rx_;
            obj.Cx_ = Cx_;            
        end
        
        
        function g = make_gauss(obj, X)
            % Instantiate the Gaussian models by providing X
            %
            %   g = obj.make_gauss(X);
            %       The distribution of Y is conditioned on X.
            %       Providing X gives rise to a set of Gaussian
            %       models of Y, whose mean vectors are A * X.
            %
            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X, 1) == r)
                error('gausslin:make_gauss:invalidarg', ...
                    'X should be an r x m real matrix.');
            end
            
            A_ = obj.A;
            
            if isempty(A_)
                AX = X;
            else
                AX = A_ * X;
            end                
                        
            if size(X,2) == 1
                mu_ = obj.mu + AX;
            else
                mu_ = bsxfun(@plus, obj.mu, AX);
            end
            
            g = gaussd_cp.from_mean_and_icov(mu_, isigma);
        end
        
        
        function P = prob(obj, X, Y)
            % Computes the likelihood on the given samples
            %
            %   P = obj.prob(X, Y);
            %       Evaluates the probability density function of the 
            %       likelihood on pairs of samples given by X and Y.
            %
            %       Here, X can be a real matrix of size r x m, and
            %       Y can be a real matrix of size d x n.            
            %       In the output, P is a real matrix of size m x n,
            %       where P(i,j) is the conditional pdf value of 
            %       Y(:,j) with respect to X(:,i).
            %
            
            P = exp(logprob(obj, X, Y));
        end        
        
        function LP = logprob(obj, X, Y)
            % Computes the logarithm-likelihood on the given samples
            %
            %   LP = obj.logprob(X, Y);
            %       Evaluates the logarithm of probability density 
            %       function of the likelihood on pairs of samples given 
            %       by X and Y.
            %
            %       Here, X can be a real matrix of size r x m, and
            %       Y can be a real matrix of size d x n.            
            %       In the output, P is a real matrix of size m x n,
            %       where P(i,j) is the conditional log-pdf value of 
            %       Y(:,j) with respect to X(:,i).
            %
                                               
            g = make_gauss(obj, X);
            LP = logprob(g, Y);
        end    
        
    end
    
    
end

