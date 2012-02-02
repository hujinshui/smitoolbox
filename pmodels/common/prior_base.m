classdef prior_base
    % The base class of prior distributions
    %
    
    % Created by Dahua Lin, on Dec 27, 2011
    %
    
    methods(Abstract)
        
        n = query_samples(obj, X);
        % Verify the validity of input samples and return the number
        
        L = logpdf(obj, X);
        % Evaluate the log pdf at given samples
        
        X = sample(obj, n);
        % Samples from the prior distribution
        
        X = pos_sample(obj, S, n);
        % Draw from posterior distribution with the stats of observations
        
        X = mapest(obj, S);
        % Performs MAP estimation with the stats of observations
        
    end    
    
end
