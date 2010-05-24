classdef gaussd_mle_pb < smodelestpb
    % The class to represent a Gaussian maximum likelihood estimation
    % problem
    %
    
    % Created by Dahua Lin, on Apr 21, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        cov_est;        % the covariance estimator
        tie_cov;       % whether to tied covariance estimation
        zero_mean;      % whether all means are zeros
        make_cp;        % whether to converted to estimated model to 
                        % gaussd_cp
    end
    
    
    methods
        
        function obj = gaussd_mle_pb(est, X)
            % Constructs a Gaussian MLE estimation problem
            %
            %   obj = gaussd_mle_pb(est, X);
            %       constructs a Gaussian MLE estimation problem.
            %
            %       In the input, est is the associated estimator of
            %       class gaussd_mle, and X is the sample matrix.
            %
            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
                error('gaussd_mle_pb:invalidarg', ...
                    'X should be a real matrix.');
            end
            
            d = est.cov_est.dim;
            if ~isempty(d)
                if size(X, 1) ~= d
                    error('gaussd_mle_pb:invalidarg', ...
                        'The dimension of X does not match that specified in the cov_est.');
                end
            end            
            
            obj = obj@smodelestpb([], X, size(X, 2));
            
            obj.cov_est = est.cov_est;
            obj.tie_cov = est.tie_cov;
            obj.zero_mean = est.zero_mean;
            obj.make_cp = est.make_cp;
        end
        
        
        function M = estimate(obj, W)
            % Perform MLE estimation of Gaussian models
            %
            %   M = obj.estimate();
            %   M = obj.estimate(W)
            %
            
            % verify input
            
            if nargin < 2
                W = [];
            end
            
            if ~isempty(W)
                if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,2) == obj.nobs)
                    error('gaussd_mle_pb:estimate:invalidarg', ...
                        'W should be either empty or a numeric real matrix with n columns.');
                end
            end
            
            X = obj.observation;
                          
            % main
            
            zmu = obj.zero_mean;
            if zmu
                mu = 0;
            else
                mu = vecmean(X, W);
            end
            
            cest = obj.cov_est;
            
            if isempty(W)
                C = cest.estcov(X, [], mu);
                m = 1;
            else                
                m = size(W, 1);
                if m == 1
                    C = cest.estcov(X, W / sum(W), mu);
                else
                    if ~obj.tie_cov
                        nw = bsxfun(@times, W, 1 ./ sum(W, 2));
                        C = cest.estcov(X, nw, mu);
                    else                        
                        wi = sum(W, 1); 
                        wi = wi / sum(wi);
                        wk = sum(W, 2);
                        wk = wk' / sum(wk);
                        C = cest.estcov_tied(X, wi, wk, mu);
                    end
                end
            end
        
            if zmu
                mu = zeros(size(X, 1), m);
            end
            
            if ~obj.make_cp
                M = gaussd_mp(mu, C);
            else
                M = gaussd_cp.from_mean_and_cov(mu, C);
            end            
        end
        
        
        function L = eval_logpri(obj, M) %#ok<INUSD,MANU>
            % Evaluate the log-prior of the estimated models
            %
            %   L = eval_logpri(obj, M);
            %
            
            L = 0;
        end
        
        
        function L = eval_loglik(obj, M) 
            % Evaluates the log-likelihood of samples 
            %
            %   L = obj.eval_loglik(M, X);
            %
            
            M = check_models(obj, M);
            L = logprob(M, X);
        end
        
        
        function M = check_models(obj, M) %#ok<MANU>
            % Check the integrity of the input models
            %
            
            if isa(M, 'gaussd_mp')
                M = to_cp(M);
            else
                if ~isa(M, 'gaussd_cp')
                    error('gaussd_mle_pb:invalidarg', ...
                        'The input models are invalid.');
                end
            end
        end
                    
    end
    
end

