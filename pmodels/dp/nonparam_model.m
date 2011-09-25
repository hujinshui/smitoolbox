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
        nobs;               % the number of observations            
        supp_inheritance;   % whether the model supports inheritance
    end
    
    %% Constructor
    
    methods
        
        function model = nonparam_model(nobs)
            % Construct the base of a nonparametric model
            %
            %   model = nonparam_model(nobs);
            %
            
            model.nobs = nobs;
        end
        
    end
    
    
    %% Interfaces
    
    methods(Abstract)
        
        a = posterior_atom(model, I, pri_a);
        % Obtain a new atom based on a specified subset of observations
        %
        %   a = model.derive_atom(I);
        %   a = model.derive_atom(I, pri_a);
        %
        %       Estimates or samples a new atom from the posterior
        %       distribution conditioned on a subset of samples, selected
        %       by the index vector I.        
        %
        %       If pri_a is omitted, the derivation is with respect to
        %       the prior given by the underlying base distribution.
        %       Otherwise, the derivation is with respect to the prior
        %       given by pri_a.
        %       
        %
        
        L = evaluate_loglik(model, a);
        % Evaluate the log-likelihood with respect to an atom
        %
        %   L = model.evaluate_loglik(a);  
        %
        %       Evaluates the log likelihood value w.r.t. the atom a.
        %
        %   L = model.evaluate_loglik(pri_a);
        %   L = model.evaluate_loglik([]);
        %
        %       Evaluates the logarithm of marginal likelihood values
        %       w.r.t. to the prior given by either pri_a, or the base
        %       distribution.               
        %
        %   L = model.evaluate_loglik(a, I);
        %   L = model.evaluate_loglik(pri_a, I);
        %   L = model.evaluate_loglik([], I);
        %
        %       Performs the evaluation for a subset of samples selected
        %       by I.
        %                
        
    end
end


