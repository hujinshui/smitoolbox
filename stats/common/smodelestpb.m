classdef smodelestpb
    % The abstract class for statistical model estimation problem
    %
    
    % Created by Dahua Lin, on Apr 19, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')        
        prior;          % the prior model of the parameters to be estimated
        observation;    % the observation        
        nobs;           % the number of observations
    end
    
    
    methods
        
        function obj = smodelestpb(pri, obs, n)
            % Constructs the estimation problem
            %
            %   obj = smodelestpb(pri, obs, n);
            %       constructs a model estimation problem that encapsulates
            %       the prior pri, and observation obs.
            %
            %       n is the number of observations.
            %        
            
            obj.prior = pri;
            obj.observation = obs;
            obj.nobs = n;
        end
    end
    
    methods(Abstract)
                
        M = estimate(obj, W);
        % Performs the estimation
        %
        %   M = obj.estimate();
        %       performs the estimation in which each observed
        %       sample is with a unit weight.
        %
        %   M = obj.estimate(W);
        %       performs the estimation with weighted samples.
        %       W can be a row vector with W(i) being the weight
        %       of the i-th sample. W can also be a matrix of
        %       size m x n, where n is the number of samples.
        %       In this case, W(k, :) gives the weights for estimating
        %       the k-th model.
        %
        %       If there are n samples, then W can be in either of
        %       the following forms:
        %       - W can be empty. In this case, each sample is with
        %         a unit weight.
        %       - W can be a row vector that gives the weights of
        %         all samples respectively.
        %       - W can be a m x n matrix where W(k, :) gives the
        %         weights for estimating the k-th model. In this case,
        %         m models are estimated, and M should be an object
        %         containing m models respectively estimated based on
        %         different sets of weights.
        %
        
        L = eval_logpri(obj, M)
        % Evaluate the log-prior of the estimated models
        %
        %   L = obj.eval_logpri(M);
        %       evaluates the log-prior of the models encapsulted in M.
        %       In particular, if obj.prior is not empty, and M is an 
        %       object that contains m models, then L will be a vector
        %       of size 1 x m, where L(i) is the log-prior of the i-th
        %       model. 
        %       If obj.prior is empty (no prior is specified), then it
        %       should simply return a scalar 0, regardless of how many
        %       models that M contains.
        %        
        
        L = eval_loglik(obj, M)
        % Evaluates the log-likelihood of samples w.r.t. estimated models
        %
        %   L = obj.eval_loglik(M);
        %       evaluates the log-likelihood of the samples
        %       with respect to the models given in M.
        %
        %       Suppose M is an object that contains m models, and there
        %       exist n samples. Then L is a real matrix of size m x n, 
        %       where L(i, j) is the log-likelihood of the j-th sample 
        %       with respect to the i-th model.
        %
                
    end
        
end