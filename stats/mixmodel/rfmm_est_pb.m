classdef rfmm_est_pb < handle
    % The class to represent a robust finite mixture model estimation problem
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        % configurational properties
        qlabeler;       % the labeler for assigning samples to component models        
        glabeler;       % the labeler for inferring inlier probabilities
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
        gprior;         % the estimated prior parameter of inlier judgment
        gposterior;     % the posterior probabilities of being inlier
        models;         % the estimated models
        Lmap;           % the likelihood map
        Lcomb;          % the combined likelihood
        Loutlier;       % the outlier likelihood
        
        modelest_pb;    % the model estimation problem
        label_pb;       % the labeling problem using current label param.
        gate_pb;        % the gating problem using current gating param.
        
        rtstatus;       % the run-time status
                
        % objective function
        
        log_mprior = nan;     % the sum of log prior value of models
        log_qprior = nan;     % the log prior of the labeling parameter        
        qscore = nan;         % the score of labeling
        log_gprior = nan;     % the log prior of the gating parameter
        gscore = nan;         % the score of gating (inlier assignment)
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
        gpri_updated;   % raised each time when gprior is updated
        gpos_updated;   % raised each time when gposterior is updated
        models_updated; % raised each time when models are updated
        
        completed;      % raised when the procedure completes
    end
    
    
    methods
        
        %% Initial setup
        
        function obj= rfmm_est_pb(est, obs, initQ, initG)
            % Constructs a finite mixture model estimation problem
            %
            %   obj = rfmm_est_pb(est, obs, initQ, initG);
            %       constructs a finite mixture model estimation problem
            %       from estimator est and observation obs.
            %
            %       Here, est is an object of class fmm_est, and 
            %       obs is the observation in a form that can be
            %       recognized by est.model_est.                        
            %
            
            % verify input arguments
                    
            if ~isa(est, 'rfmm_est')
                error('rfmm_est_pb:invalidarg', ...
                    'est should be an object of est.');
            end
            
            if ~(isfloat(initQ) && isreal(initQ) && ndims(initQ) == 2)
                error('rfmm_est_pb:invalidarg', ...
                    'initQ should be a real matrix.');
            end
               
            if ~(isfloat(initG) && isreal(initG) && ndims(initG) == 2 && size(initG,1) == 1)
                error('rfmm_est_pb:invalidarg', ...
                    'initG should be a real scalar or real row vector.');
            end
            
            [m, n] = size(initQ);
                        
            % set fields
            
            obj.qlabeler = est.qlabeler;
            obj.glabeler = est.glabeler;
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
                error('rfmm_est_pb:invalidarg', ...
                    'The number of observations does not match the size of initQ.');
            end
            
            if isscalar(initG)
                obj.gposterior = initG(ones(1, n));
            else
                if size(initG, 2) ~= n
                    error('rfmm_est_pb:invalidarg', ...
                        'The size of initG is not consistent with the number of samples.');
                end
                obj.gposterior = initG;
            end
            
            obj.Loutlier = est.lp_outlier(ones(1, n));
           
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
                                
                update_gpri(obj);
                
                if elvl >= 3
                    notify(obj, 'gpri_updated');
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
                
                update_gpos(obj);
                
                if elvl >= 3
                    notify(obj, 'gpos_updated');
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
        
        function update_gpri(obj)
            % Updates the prior parameter of inlier assignment
            
            gl = obj.glabeler;
            G = obj.gposterior;
            
            gpri = gl.estimate_param(G);
            gpb = gl.accept_param(gpri);
            
            obj.gprior = gpri;
            obj.gate_pb = gpb;
            
            obj.log_gprior = gl.eval_logpri(gpri);
            obj.gscore = gpb.eval_qscore(G);
            
            update_total_objv(obj);
        end
                        
        
        function update_qpos(obj)
            % Update the posterior distribution of labels
            
            qpb = obj.label_pb;
            L = obj.Lmap;
            gL = bsxfun(@times, obj.gposterior, L);
            Q = qpb.infer_q(gL);
            obj.qposterior = Q;
            
            obj.qscore = qpb.eval_qscore(Q);
            obj.Lcomb = sum(Q .* L, 1);
            
            update_sumloglik(obj);            
            update_total_objv(obj);
        end
        
        
        function update_gpos(obj)
            % Update the posterior probabilities of inlier
            
            gpb = obj.gate_pb;
            
            L1 = obj.Lcomb;
            L0 = obj.Loutlier;
            G = gpb.infer_q(L1, L0);
            obj.gposterior = G;
            
            obj.gscore = gpb.eval_qscore(G);
            
            update_sumloglik(obj);
            update_total_objv(obj);            
        end
                
        
        function update_models(obj)
            % Update the observation model parameters
            
            mpb = obj.modelest_pb;
            Q = obj.qposterior;            
            W = bsxfun(@times, obj.gposterior, Q);
            mdls = mpb.estimate(W);
            
            obj.models = mdls;
            L = mpb.eval_loglik(mdls);
            
            obj.Lmap = L;
            obj.Lcomb = sum(Q .* L, 1);
            
            obj.log_mprior = sum(mpb.eval_logpri(mdls));
            
            update_sumloglik(obj);            
            update_total_objv(obj);
        end
                             
    end
    
    
    methods(Access='private')
        
        function update_sumloglik(obj)
            
            g = obj.gposterior;            
            L1 = obj.Lcomb;
            L0 = obj.Loutlier;
            
            v = sum(g .* L1 + (1-g) .* L0);
            
            obj.sum_loglik = v;
        end        
        
        function update_total_objv(obj)
            
               obj.total_objv = ...
                obj.log_mprior + ...
                obj.log_qprior + ...
                obj.log_gprior + ...
                obj.qscore + ...
                obj.gscore + ...
                obj.sum_loglik;     
        end
        
    end

end

