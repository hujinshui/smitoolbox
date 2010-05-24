classdef multd_est_pb < smodelestpb
    % The class to represent a multi-class discrete distribution estimation 
    % problem
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
        
    properties
        dim;        % the dimension of the p-space (# classes)
    end
    
    
    methods        
        function obj = multd_est_pb(est, Q)
            % Constructs a binary distribution estimation problem
            %
            %   obj = multd_map(est, Q);
            %       constructs a maximum likelihood estimation problem
            %       
            %       In the input, est is the associated estimator, and
            %       Q is the soft label map for samples.
            %
                                    
            if ~(isfloat(Q) && isreal(Q) && ndims(Q) == 2)
                error('multd_est_pb:invalidarg', ...
                    'Q should be a real matrix.');
            end     
            
            d = size(Q, 1);
            pri = est.prior;            
            if ~isempty(pri)
                if d ~= pri.dim
                    error('multd_est_pb:invalidarg', ...
                        'The dimension of Q does not match that of prior.');
                end
            end                            
            
            obj = obj@smodelestpb(pri, Q, size(Q, 2));
            obj.dim = d;
        end
        
        
        function p = estimate(obj, W) 
            % Estimate the Gaussian mean in MAP manner
            %
            %   p = obj.estimate();
            %   p = obj.estimate(W);
            %       estimates the binary distribution p from the 
            %       soft assignment of each sample.            
            %
            %       In the output, p is matrix of size d x m, where
            %       m is the number of models to estimate
            %
            
            % verify input
                        
            n = obj.nobs;
            
            if nargin >= 2 && ~isempty(W)
                if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,2) == n)
                    error('multd_est_pb:estimate:invalidarg', ...
                        'W should be either empty or a real matrix with n columns.');
                end
            end
            
            % main
            
            Q = obj.observation;
           
            if isempty(W)
                q = sum(Q, 2);
            else
                q = Q * W';
            end
            
            pri = obj.prior;
            if ~isempty(pri)
                q = bsxfun(@plus, q, pri.alpha - 1);
            end
            
            p = bsxfun(@times, q, 1 ./ sum(q, 1));     
        end
        
        
        function lpri = eval_logpri(obj, p)
            % Evaluates the log-prior of a binary distribution
            %
            %   lpri = eval_logpri(obj, p);
            %                        
            
            pri = obj.prior;
            if isempty(pri)
                lpri = 0;
            else
                lpri = pri.logprob(p);
            end
        end
        
        
        function L = eval_loglik(obj, p)
            % Evaluates the log-likelihood of observed samples
            %
            %   L = eval_loglik(obj, p);
            %
            
            d = obj.dim;
            
            if ~(isfloat(p) && isreal(p) && ndims(p) == 2 && size(p, 1) == d)
                error('multd_est_pb:eval_loglik:invalidarg', ...
                    'p should be a real matrix with size(p,1) == d.');
            end
            
            Q = obj.observation;            
            lp = log(p);
            lp(p == 0) = 0;
            L = lp' * Q;
        end
        
    end    
end

