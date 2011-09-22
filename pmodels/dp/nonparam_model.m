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
        
        a = create_atom(model, I, pri_a, t);
        % Obtain a new atom based on a specified subset of observations
        %
        %   a = model.estimate_atom(I);
        %
        %       estimates/samples a new atom conditioned on a subset of
        %       observations, and with respect to the base prior. 
        %
        %       Here, I is a vector of indices. If I is empty, it samples
        %       from the prior distribution.
        %
        %   For the model that supports inheritance:
        %
        %   a = model.estimate_atom(I, pri_a, t);
        %
        %       estimates/samples a new atom conditioned on a subset of
        %       observations, and with respect to the prior formed via 
        %       t-fold transition of a prior atom a.
        %
        
        L = evaluate_loglik(model, a, t);
        % Evaluate the log-likelihood with respect to an atom
        %
        %   L = model.evaluate_loglik(a);
        %
        %       evaluates the log-likelihood of all observations with
        %       respect to the input atom a.
        %
        %       In the output, L is an 1 x nobs row vector.
        %
        %   L = model.evaluate_loglik();
        %
        %       evaluates the logarithm of marginal likelihood of all 
        %       observations, with respect to the base prior.                
        %
        %   For the model that supports inheritance:
        %
        %   L = model.evaluate_loglik(a, t);
        %
        %       evaluates the logarithm of marginal likelihood of all
        %       observations, with respect to the prior formed via
        %       t-fold transition of the atom a.
        %
        
        lpri = evaluate_logpri(model, a, pri_a, t);
        % Evaluates the log-prior of an atom
        %
        %   lpri = model.evaluate_logpri(a);
        %       
        %       evaluates the log-prior of the input atom a, with 
        %       respect to the base prior.
        %
        %   For the model that supports inheritance
        %
        %   lpri = model.evaluate_logpri(a, pri_a, t);
        %
        %       evaluates the log-prior of the input atom a, with
        %       respect to the prior formed by t-fold transition of
        %       the prior atom (pri_a).
        %
        
        
    end
end
