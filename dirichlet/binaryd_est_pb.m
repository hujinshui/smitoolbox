classdef binaryd_est_pb < smodelestpb
    % The class to represent a binary distribution estimation problem
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
        
    
    methods        
        function obj = binaryd_est_pb(est, Q)
            % Constructs a binary distribution estimation problem
            %
            %   obj = binaryd_map(est, Q);
            %       constructs a maximum likelihood estimation problem
            %       
            %       In the input, est is the associated estimator, and
            %       Q is the vector of soft assignment to label 1.
            %
            
            if ~(isfloat(Q) && isreal(Q) && ndims(Q) == 2 && size(Q,1) == 1)
                error('binaryd_est_pb:invalidarg', ...
                    'Q should be a real row vector.');
            end     
            
            obj = obj@smodelestpb(est.prior, Q, size(Q, 2));
        end
        
        
        function p = estimate(obj, W) 
            % Estimate the Gaussian mean in MAP manner
            %
            %   p = obj.estimate();
            %   p = obj.estimate(W);
            %       estimates the binary distribution p from the 
            %       soft assignment of each sample.            
            %
            %       In the output, p is a scalar or a 1 x m row vector
            %       if there are m rows in W.
            %
            
            % verify input
                        
            n = obj.nobs;
            
            if nargin >= 2 && ~isempty(W)
                if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,2) == n)
                    error('binaryd_est_pb:estimate:invalidarg', ...
                        'W should be either empty or a real matrix with n columns.');
                end
            end
            
            % main
            
            Q = obj.observation;
           
            if isempty(W)
                q1 = sum(Q, 2);
                q0 = n - q1;
            else
                q1 = Q * W';
                q0 = sum(W, 2).';
            end
            
            pri = obj.prior;
            if ~isempty(pri)
                q1 = q1 + (pri.alpha - 1);
                q0 = q0 + (pri.beta - 1);
            end
            
            p = q1 ./ (q1 + q0);      
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
            
            if ~(isfloat(p) && isreal(p) && ndims(p) == 2 && size(p, 1) == 1)
                error('binaryd_est_pb:eval_loglik:invalidarg', ...
                    'p should be a real row vector.');
            end
            
            lp1 = log(p).';
            lp1(p == 0) = 0;
            
            lp0 = log(1 - p).';
            lp0(p == 1) = 0;
            
            Q = obj.observation;            
            L = lp1 * Q + lp0 * (1 - Q);
        end
        
    end    
end

