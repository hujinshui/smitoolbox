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
        
        Zmethod;    % the method code for choosing Z
                    % 0 - let Z = Q, the posterior probabilities
                    % 1 - set Z to the best label
                    % 2 - sample Z
    end
    
    % observations
    
    properties(GetAccess='public', SetAccess='private')
        obs;        % the observations
        nobs;       % the number of observations (n)
    end
    
    
    % run-time state
    
    properties
        K;          % the number of mixture components (K)
        Pi;         % the distribution over mixture components [K x 1]
        params;     % the estimated component parameters
        Z;          % the posterior probabilities or the estimated labels
        
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
            %                   - 'hard-em':    hard E-M method, assigning
            %                                   each sample completely
            %                                   to the best class
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
                    zmd = 0;
                    samp = false;
                case 'hard-em'
                    zmd = 1;
                    samp = false;
                case 'gibbs'
                    zmd = 2;
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
            obj.Zmethod = zmd;
            obj.sampling = samp;
        end
        
    end
    
    %% Interface methods
    
    methods
        
        function obj = initialize(obj, obs, params_init)
            % Initialize the FMM estimator state
            %
            %   obj = obj.initialize(obs, params_init);
            %       initializes the state of finite mixture model
            %       estimator.
            %
            %       Input arguments:
            %       - obs:          the observations
            %       - params_init:  the initial set of component parameters
            %
            
            gm = obj.gmodel;
            
            n = gm.query_obs(obs);
            K_ = gm.query_params(params_init);
            
            obj.obs = obs;
            obj.nobs = n;
            
            obj.K = K_;
            obj.Pi = (1/K_) * ones(K_,1);
            obj.params = params_init;            
            obj.Llik = gm.loglik(params_init, obs);
        end
        
        
        function obj = update(obj)
            % Update the state
            %
            %   obj.update();
            %       updates the object state (perform one-step of E-M
            %       optimization or one move of Gibbs sampling)
            %
            
            gm = obj.gmodel;
            K_ = obj.K;
            
            % E-step (labeling)
            
            Pi_ = obj.Pi;
            L = obj.Llik;
            
            Q = ddposterior(Pi_, L, 'LL');
            
            zmd = obj.Zmethod;
            if zmd == 0
                Z_ = Q;
            elseif zmd == 1
                [~, Z_] = max(Q, [], 1);
            else
                Z_ = ddsample(Q, 1);
            end
            
            obj.Z = Z_;
            
            % M-step (parameter estimation)
            
            % estimate component parameters
            
            pri = obj.prior;
            X = obj.obs;
            samp = obj.sampling;
            
            if zmd == 0
                V = Z_;
            else
                V = {K_, Z_};
            end
            
            if isempty(pri)
                thetas = gm.mle(X);
            else
                cap = gm.capture(X, V);
                
                if samp
                    thetas = pri.pos_sample(cap);
                else
                    thetas = pri.mapest(cap);
                end
            end
            
            obj.params = thetas;
            obj.Llik = gm.loglik(thetas, X);
            
            % estimate component prior
            
            if samp                
                H = accum_counts(V);
                obj.Pi = dird_pos_sample(obj.pricount + 1, H);
            else
                obj.Pi = ddestimate(V, [], obj.pricount);
            end
            
        end
        
        
        function R = output(obj)
            % Outputs a sample
            %
            %   R = obj.output();
            %
            
            R.K = obj.K;
            R.Pi = obj.Pi;
            R.params = obj.params;
            R.Z = obj.Z;
            
        end
        
        
        function b = is_ready(obj)
            % Tests whether the object is ready for running
            %
            %   b = obj.is_ready();
            %
            
            b = ~isempty(obj.obs);            
        end
        
    end
    
    
    %% Objective evaluation
    
    methods
        
        function objv = evaluate_objv(obj)
            % Evaluate the objective function of the current state
            %
            %   objv = obj.evaluate_objv();
            %
            
            % log-pri: Pi
            
            c0 = obj.pricount;
            
            log_pi = log(obj.Pi);
            
            if isequal(c0, 0)
                lpri_pi = 0;
            else
                lpri_pi = c0 * sum(log_pi);
            end
            
            % log-pri: params
            
            thetas = obj.params;            
            pri = obj.prior;
            
            if isempty(pri)
                lpri_t = 0;
            else
                lpri_t = pri.logpdf(thetas);
                lpri_t = sum(lpri_t);
            end
            
            % log-lik: labels (Z)
            
            zmd = obj.Zmethod;
            Z_ = obj.Z;
            
            if zmd == 0
                llik_z = sum(log_pi' * Z_);
            else
                llik_z = sum(log_pi(Z_));
            end
            
            % log-lik: observations
            
            L = obj.Llik;
            
            if zmd == 0
                llik_x = sum(Z_ .* L, 1);
            else                
                n = size(L, 2);
                llik_x = L(sub2ind(size(L), 1:n, Z_));
            end
            llik_x = sum(llik_x);
            
            % entropy
            
            if zmd == 0
                ent = ddentropy(Z_);
            else
                ent = 0;
            end
                
            % overall sum
            
            objv = lpri_pi + lpri_t + llik_z + llik_x + ent;
            
        end        
    
    end
    
    
    
end




