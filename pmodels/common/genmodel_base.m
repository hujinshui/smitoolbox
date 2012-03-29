classdef genmodel_base
    % The base of the classes that implement simple generative model
    %
    %   This is an abstract class that defines the functions to be
    %   implemented by derived classes.
    %
    
    methods(Abstract)         
        
        n = query_obs(obj, X);
        % Verify the observations and return the number of observations
        
        n = query_params(obj, A);
        % Verify the parameters and return the number of parameters
                
        L = loglik(obj, A, X);
        % Evaluate the log-likelihood of all samples w.r.t all params
        
        A = mle(obj, X, W, I);
        % Performs maximum likelihood estimation of parameters
        %
        % W: empty or a matrix of size n x K
        %
                
        S = capture(obj, X, W, I);
        % Captures the sufficient stats of observations as updates to prior
        
    end
    
end
