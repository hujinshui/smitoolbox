classdef nonparam_model
    % The base class of all non-parametric models
    %
    %   This class is an abstract class that defines the interfaces
    %   for a non-parametric model. 
    %
    
    % Created by Dahua Lin, on Sep 17, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='protected')                  
        supp_inheritance = false;   % whether the model supports inheritance
    end    
    
    
    %% Interfaces
    
    methods(Abstract)
        
        n = get_num_samples(model, X);
        % Verifies the validity of a data-set, and returns #samples
        %
        %   a = model.get_num_samples(X);
        %
        %       Verifies the validity of X as a dataset, and returns the
        %       number of samples in X.
        %
        
        
        a = posterior_atom(model, X, I, pri_a);
        % Obtain a new atom based on a specified subset of observations
        %
        %   a = model.posterior_atom(X, I);
        %   a = model.posterior_atom(X, I, pri_a);
        %
        %       Estimates or samples a new atom from the posterior
        %       distribution conditioned on a subset of samples in X, 
        %       selected by the index vector I.        
        %
        %       If pri_a is omitted, the derivation is with respect to
        %       the prior given by the underlying base distribution.
        %       Otherwise, the derivation is with respect to the prior
        %       given by pri_a.
        %       
        %
        
        L = evaluate_loglik(model, a, X);
        % Evaluate the log-likelihood with respect to an atom
        %
        %   L = model.evaluate_loglik(a, X);  
        %
        %       Evaluates the log likelihood value w.r.t. the atom a.
        %
        %   L = model.evaluate_loglik(pri_a, X);
        %   L = model.evaluate_loglik([], X);
        %
        %       Evaluates the logarithm of marginal likelihood values
        %       w.r.t. to the prior given by either pri_a, or the base
        %       distribution.                       
        %                
        
    end
end


