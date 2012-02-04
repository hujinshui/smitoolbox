classdef fmm_std < smi_state
    % The class that implements standard finite mixture model
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Sep 27, 2011
    %       - Modified by Dahua Lin, on Dec 26, 2011
    %
    
    %% Properties
    
    % configurations
    
    properties(GetAccess='public', SetAccess='private')
        gmodel;     % The underlying generative model
        prior;      % The parameter prior
        pricount;   % The prior count of each component                                    
        sampling;   % whether it uses sampling
    end
    
    % observations
    
    properties(GetAccess='public', SetAccess='private')
        obs;        % the observations
        nobs;       % the number of observations (n)
        weights;    % the weights of observations (empty or n x 1)
    end
    
    
    % run-time state
    
    properties
        K;          % the number of mixture components
        sol;        % the current solution (a struct with fields)                
        Llik;       % the log-likelihood w.r.t. all components [K x n]
    end
    
    
    %% Construction
    
    methods
        
        function obj = fmm_std(method, gm, pri, c0)
            % Create a standard finite mixture model estimator
            %
            %   obj = fmm_std(method, gm, pri);
            %   obj = fmm_std(method, gm, pri, c0);
            %
            %       Creates a standard finite mixture model estimator
            %
            %       Inputs:
            %       - method:   either of the following method names:
            %                   - 'em':         standard E-M method
            %                   - 'gibbs':      Gibbs sampling
            %
            %       - gm:       the generative model object
            %
            %       - pri:      the prior model object or empty
            %
            %       - c0:       the prior count of each component
            %                   (if omitted, it is set to zero).
            %
            
            % verify input
            
            if ~ischar(method)
                error('fmm_std:invalidarg', 'method should be a char string.');
            end
            
            switch lower(method)
                case 'em'
                    samp = false;
                case 'gibbs'
                    samp = true;
                otherwise
                    error('fmm_std:invalidarg', ...
                        'Invalid method name %s', method);
            end
            
            if ~isa(gm, 'genmodel_base')
                error('fmm_std:invalidarg', ...
                    'gm should be an instance of a sub-class of genmodel_base.');
            end
            
            if isempty(pri)
                if samp
                    error('fmm_std:invalidarg', ...
                        'pri should be provided when using sampling.');
                end
            else
                if ~isa(pri, 'prior_base')
                    error('fmm_std:invalidarg', ...
                        'pri should be an instance of a sub-class of prior_base.');
                end
            end
            
            
            if nargin < 4
                c0 = 0;
            else
                if ~(isfloat(c0) && isreal(c0) && isscalar(c0) && c0 >= 0)
                    error('fmm_std:invalidarg', ...
                        'c0 should be a nonnegative real value scalar.');
                end
            end
            
            % set fields
            
            obj.gmodel = gm;
            obj.prior = pri;
            obj.pricount = double(c0);
            obj.sampling = samp;
        end
        
    end
    
    %% Interface methods
    
    methods
        
        function obj = initialize(obj, X, w, method, arg)
            % Initialize the FMM estimator state
            %
            %   obj = obj.initialize(X, w, 'params', A);
            %
            %       initialize with given initial parameters
            %
            %   obj = obj.initialize(X, w, 'labels', z);
            %
            %       initialize with given initial labels
            %
            %   obj = obj.initialize(X, w, 'Q', Q);
            %   
            %       initialize with soft assignment matrix
            %
            %   obj = obj.initialize(X, w, 'rand', K);
            %
            %       initialize randomly with K component mixtures.
            %
            %   Here, X is a sample matrix of size d x n, and 
            %   w is either empty (all samples have a unit weight), or
            %   a vector of length n.
            %
            
            gm = obj.gmodel;
            s = fmm_init(obj.prior, gm, X, w, method, arg);
           
            if obj.sampling
                s.z = [];
            else
                s.Q = [];
            end
            
            obj.obs = X;
            obj.nobs = gm.query_obs(X);
            if ~isempty(w)
                if size(w, 2) > 1; w = w.'; end
                obj.weights = w;
            end
            
            obj.K = s.K;
            obj.sol = s;       
            obj.Llik = gm.loglik(s.params, X);
        end        
               
        
        function obj = update(obj)
            % Update the state
            %
            %   obj.update();
            %       updates the object state (perform one-step of E-M
            %       optimization or one move of Gibbs sampling)
            %
            
            pri = obj.prior;
            gm = obj.gmodel;
            X = obj.obs;
            w = obj.weights;
            c0 = obj.pricount;
            
            if obj.sampling                
                [obj.sol, obj.Llik] = fmm_gibbs_update( ...
                    pri, gm, X, w, c0, obj.sol, obj.Llik);                
            else
                [obj.sol, obj.Llik] = fmm_em_update(...
                    pri, gm, X, w, c0, obj.sol, obj.Llik);
            end
            
        end
        
        
        function s = output(obj)
            % Outputs a sample
            %
            %   s = obj.output();
            %
            
            s = obj.sol;            
        end
        
        
        function b = is_ready(obj)
            % Tests whether the object is ready for running
            %
            %   b = obj.is_ready();
            %
            
            b = ~isempty(obj.sol);            
        end
        
        
        function s = merge_samples(obj, samples)
            % Merges multiple samples into an optimal one
            %
            %   s = obj.merge_samples(samples);
            %
            
            s = fmm_merge_samples(obj.prior, obj.gmodel, obj.pricount, ...
                obj.obs, obj.weights, samples);            
        end        
        
    end
    
    
    %% Objective evaluation
    
    methods
        
        function objv = evaluate_objv(obj)
            % Evaluate the objective function of the current state
            %
            %   objv = obj.evaluate_objv();
            %
            
            if ~obj.sampling
                objv = fmm_em_objective(...
                    obj.prior, obj.Llik, obj.weights, obj.pricount, obj.sol);
            else
                error('fmm_std:invalidarg', ...
                    'evaluate_objv can only be invoked in non-sampling mode.');
            end
            
        end        
    
    end          
    
    
end




