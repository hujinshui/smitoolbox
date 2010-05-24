classdef gaussd_map_pb < smodelestpb
    % The class to represent a Gaussian Maxium-a-Posteriori problem
    %
    
    % Created by Dahua Lin, on April 21, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        icov;       % the inverse covariance of observations
    end
    
    methods
        
        function obj = gaussd_map_pb(est, X)
            % Constructs an Gaussian MAP estimation problem
            %
            %   obj = gaussd_map_pb(est, X);
            %       In the input, est is the associated estimator of
            %       class gaussd_mle, and X is the sample matrix.
            %
            
            d = est.dim;
            
            assert(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == d, ...
                'gaussd_map_pb:invalidarg', ...
                'X should be a real matrix of size d x n.');
            
            obj = obj@smodelestpb(est.prior, X, size(X, 2));            
            obj.icov = est.icov;
        end
        
        
        function mu = estimate(obj, W)
            % Perform MLE estimation of Gaussian models
            %
            %   mu = obj.estimate();
            %   mu = obj.estimate(W);
            %
            
            % verify input
            
            if nargin < 2
                W = [];
            end
            
            if ~isempty(W)
                if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,2) == obj.nobs)
                    error('gaussd_map_pb:estimate:invalidarg', ...
                        'W should be either empty or a numeric real matrix with n columns.');
                end
                
                tw = sum(W, 2).';
                K = size(W, 1);
            else
                tw = obj.nobs;
                K = 1;
            end
            
            X = obj.observation;
                          
            % main
            
            isigma = obj.icov;
            g0 = obj.prior;
            
            % For theta2
            U2 = tw * isigma;
            
            % For theta1
            LX = isigma * X;
            if isempty(W)
                U1 = sum(LX, 2);
            else
                U1 = LX * W.';
            end
            
            % compute mean
            if isempty(g0)
                T1 = U1;
            else
                if K == 1
                    T1 = g0.theta1 + U1;
                else
                    T1 = bsxfun(@plus, g0.theta1, U1);
                end
            end
            
            if isempty(g0)
                T2 = U2;
            else
                T2 = g0.theta2 + U2;
            end
            
            mu = cldiv(T2, T1);                 
        end
        
        
        function L = eval_logpri(obj, mu)
            % Evaluate the log-prior of the estimated means
            %
            %   L = eval_logpri(obj, mu);
            %
            
            g0 = obj.prior;
            check_mu(obj, mu);
            
            if isempty(g0)
                L = 0;
            else
                L = logprob(g0, mu);
            end            
        end
        
        
        function L = eval_loglik(obj, mu) 
            % Evaluates the log-likelihood of samples w.r.t. estimated
            % means
            %
            %   L = obj.eval_loglik(M, X);
            %
                        
            check_mu(obj, mu);
            g = gaussd_cp.from_mean_and_icov(mu, obj.icov);
            X = obj.observation;
            L = logprob(g, X);
        end
        
        
        function check_mu(obj, mu)
            % Check whether the input mu is correct
            %
            
            d = obj.icov.dim;
            if ~(isfloat(mu) && isreal(mu) && ndims(mu) == 2 && size(mu,1) == d)
                error('gaussd_map_pb:invalidarg', ...
                    'mu should be a real matrix with size(mu,1) == d.');
            end
        end
    end
end

