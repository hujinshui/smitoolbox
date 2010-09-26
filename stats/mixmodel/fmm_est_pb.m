classdef fmm_est_pb < handle
    % The class to represent a finite mixture model estimation problem
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        % configurational properties
        qlabeler;       % the labeler for assigning samples to component models        
        model_est;      % the estimator of component models from observations
        
        observation;    % the observed samples
        nobs;           % the number of observed samples (n)
        ncomps;         % the number of component models (m)    
        
        maxiter;        % the maximum number of iterations
        tolerance;      % the tolerance of change at convergence
        event_level;    % the level of event notification
        
        % solution        
        qprior;         % the estimated prior parameter of labeling
        qposterior;     % the posterior label probabilities
        models;         % the estimated models
        Lmap;           % the likelihood map
        
        modelest_pb;    % the model estimation problem
        label_pb;       % the labeling problem using current label param.
        
        rtstatus;       % the run-time status
                
        % objective function
        
        log_mprior = nan;     % the sum of log prior value of models
        log_qprior = nan;     % the log prior of the labeling parameter
        qscore = nan;         % the score of labeling
        sum_loglik = nan;     % the sum of log-likelihood, defined by
                              % \sum_i \sum_k q(k,i) \log p(x_i|model_k) 
                        
        total_objv = nan;     % total objective value
    end
    
    
    events(NotifyAccess='private')
        % Note that these events will only be raised in the solve method
        
        iter_start;     % raised when an iteration starts
        iter_end;       % raised when an iteration ends
        
        qpri_updated;   % raised each time when qprior is updated
        qpos_updated;   % raised each time when qposterior is updated
        models_updated; % raised each time when models are updated
        
        completed;      % raised when the procedure completes
    end
    
    
    methods
        
        %% Initial setup
        
        function obj= fmm_est_pb(est, obs, initQ)
            % Constructs a finite mixture model estimation problem
            %
            %   obj = fmm_est_pb(est, obs, initQ);
            %       constructs a finite mixture model estimation problem
            %       from estimator est and observation obs.
            %
            %       Here, est is an object of class fmm_est, and 
            %       obs is the observation in a form that can be
            %       recognized by est.model_est.
            %
            
            % verify input arguments
                    
            if ~isa(est, 'fmm_est')
                error('fmm_est_pb:invalidarg', ...
                    'est should be an object of est.');
            end
            
            if ~(isfloat(initQ) && isreal(initQ) && ndims(initQ) == 2)
                error('fmm_est_pb:invalidarg', ...
                    'initQ should be a real matrix.');
            end
            
            [m, n] = size(initQ);
                        
            % set fields
            
            obj.qlabeler = est.qlabeler;
            obj.model_est = est.model_est;
            
            obj.observation = obs;                                                
            obj.nobs = n;
            obj.ncomps = m;
            
            obj.maxiter = est.maxiter;
            obj.tolerance = est.tolerance;
            obj.event_level = est.event_level;
            
            % initialize solution
            
            obj.qposterior = initQ;
            
            mpb = est.model_est.accept(obs);
            if mpb.nobs ~= n
                error('fmm_est_pb:invalidarg', ...
                    'The number of observations does not match the size of initQ.');
            end
            
            obj.modelest_pb = mpb;            
        end
        
        
        %% Main skeleton
        
        function info = solve(obj)
            % Solve the estimation problem
            %
            %   obj.solve();
            %       solves the estimation problem. The estimated results
            %       can be directly read from the object property.
            %
            %   info = obj.solve();
            %       solves the problem and additionally returns the 
            %       information relevant to the procedure.
            %
            %       info is a struct that contains the following fields:
            %       - niters:       the number of iterations elapsed
            %       - converged:    whether the procedure converges
            %       - objvs:        the curve of objective value
            %       
            
            T = obj.maxiter;            
            objvs = zeros(1, T);
            
            tol = obj.tolerance;
            elvl = obj.event_level;
                
            obj.rtstatus = struct('t', 0, 'converged', false);
            t = 0;
            
            while ~obj.rtstatus.converged && t <= T
                t = t + 1;
                obj.rtstatus.t = t;
                
                if elvl >= 2
                    notify(obj, 'iter_start');
                end
                
                % M-steps
                                                
                update_qpri(obj);
                
                if elvl >= 3
                    notify(obj, 'qpri_updated');
                end                
                
                update_models(obj);
                
                if elvl >= 3
                    notify(obj, 'models_updated');
                end
                
                % E-steps
                                
                update_qpos(obj);
                
                if elvl >= 3
                    notify(obj, 'qpos_updated');
                end 
                                
                % decide convergence
                
                objvs(t) = obj.total_objv;
                if t > 1
                    objv_ch = objvs(t) - objvs(t-1);
                    if abs(objv_ch) <= tol
                        obj.rtstatus.converged = true;
                    end
                end
                
                if elvl >= 2
                    notify(obj, 'iter_end');
                end
            end
                      
            if elvl >= 1
                notify(obj, 'completed');
            end
            
            obj.rtstatus = [];
            
            % output info
                                            
            if nargout >= 1
                info = struct( ...
                    'niters', t, ...
                    'converged', converged, ...
                    'objvs', objvs(1:t));
            end
        end
        
        %% Updating Steps
        
        function update_qpri(obj)
            % Updates the prior parameter of labels
            
            ql = obj.qlabeler;
            Q = obj.qposterior;
            
            qpri = ql.estimate_param(Q);
            qpb = ql.accept_param(qpri);
            
            obj.qprior = qpri;
            obj.label_pb = qpb;
            
            obj.log_qprior = ql.eval_logpri(qpri);
            obj.qscore = qpb.eval_qscore(Q);
            
            update_total_objv(obj);
        end
        
        function update_qpos(obj)
            % Update the posterior distribution of labels
            
            qpb = obj.label_pb;
            Q = qpb.infer_q(obj.Lmap);
            obj.qposterior = Q;
            
            obj.qscore = qpb.eval_qscore(Q);
            L = obj.Lmap;
            obj.sum_loglik = sum(sum(Q .* L, 1));
            
            update_total_objv(obj);
        end
        
        function update_models(obj)
            % Update the observation model parameters
            
            mpb = obj.modelest_pb;
            mdls = mpb.estimate(obj.qposterior);
            
            obj.models = mdls;
            L = mpb.eval_loglik(mdls);
            obj.Lmap = L;
            
            obj.log_mprior = sum(mpb.eval_logpri(mdls));
            Q = obj.qposterior;
            obj.sum_loglik = sum(sum(Q .* L, 1));
            
            update_total_objv(obj);
        end
                             
    end
    
    
    methods(Access='private')
        
        function update_total_objv(obj)
                        
            obj.total_objv = ...
                obj.log_mprior + ...
                obj.log_qprior + ...
                obj.qscore + ...
                obj.sum_loglik;
        end
        
    end

end

