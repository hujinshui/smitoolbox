classdef fmm_est_mon < handle
    % The class implementing the standard monitor of finite mixture model
    % estimation
    %
    
    % Created by Dahua Lin, on Apr 22, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        output_fid = 1;
        problem;
        event_level;
        
        % snapshot of the problem status
        log_mprior;
        log_qprior;
        qscore;
        sum_loglik;
        total_objv_start;
        total_objv;        
    end    
    
    methods
        
        function obj = fmm_est_mon(problem, eventlevel)
            % Constructs a monitor connected to a fmm estimation problem
            %
            %   obj = fmm_est_mon(problem);
            %   obj = fmm_est_mon(problem, eventlevel);
            %
            %       By default, it uses problem's event_level.             
            %
            
            assert(isa(problem, 'fmm_est_pb'), ...
                'fmm_est_mon:invalidarg', ...
                'The first argument should be an object of fmm_est_pb.');
            
            if nargin < 2
                eventlevel = problem.event_level;
            else            
                assert(isnumeric(eventlevel) && isscalar(eventlevel), ...
                    'fmm_est_mon:invalidarg', ...
                    'The eventlevel should be a numeric scalar.');
            end
            
            obj.problem = problem;
            obj.event_level = eventlevel;
            
            % connect to the problem                    

            if eventlevel >= 1
                addlistener(problem, 'completed', @(src, evnt) on_completed(obj, src));
            end
            
            if eventlevel >= 2
                addlistener(problem, 'iter_start', @(src, evnt) on_iter_start(obj, src));
                addlistener(problem, 'iter_end', @(src, evnt) on_iter_end(obj, src));
            end
            
            if eventlevel >= 3
                addlistener(problem, 'qpri_updated', @(src, evnt) on_qpri_updated(obj, src));
                addlistener(problem, 'qpos_updated', @(src, evnt) on_qpos_updated(obj, src));
                addlistener(problem, 'models_updated', @(src, evnt) on_models_updated(obj, src));
            end            
        end
    end
                        
    methods
    
        function on_iter_start(obj, problem) 
                   
            elvl = obj.event_level; 
            if elvl >= 2
                obj.total_objv_start = problem.total_objv;
                
                if elvl >= 3                    
                    fprintf('Iter %d:\n', problem.rtstatus.t);
                    obj.take_snapshot(problem);                    
                end                
            end
        end
        
        function on_iter_end(obj, problem) 
            
            elvl = obj.event_level;
            if elvl >=2                 
                t = problem.rtstatus.t;
                
                sv = obj.total_objv_start;
                v = problem.total_objv;
                obj.total_objv = v;
                
                if elvl < 3
                    fprintf(obj.output_fid, ...
                        'Iter %d: objv = %g [ch = %g]\n', t, v, v - sv);
                else
                    fprintf(obj.output_fid, ...
                        '\tSum up: objv = %g [ch = %g]\n', v, v - sv);
                end
            end
        end
        
        function on_qpri_updated(obj, problem) 
            
            if obj.event_level >= 3                
                qp0 = obj.log_qprior;
                qs0 = obj.qscore;
                
                qp = problem.log_qprior;
                qs = problem.qscore;
                
                fprintf(obj.output_fid, ...
                    '\tqprior updated: log_qprior %g --> %g, qscore %g --> %g [ch = %g]\n', ...
                    qp0, qp, qs0, qs, (qp + qs) - (qp0 + qs0));    
                
                obj.take_snapshot(problem);
            end                                                    
        end
        
        function on_qpos_updated(obj, problem)  
            
            if obj.event_level >= 3
                qs0 = obj.qscore;
                sl0 = obj.sum_loglik;
                
                qs = problem.qscore;
                sl = problem.sum_loglik;
                
                fprintf(obj.output_fid, ...
                    '\tqposterior updated: qscore %g --> %g, sum_loglik %g --> %g [ch = %g]\n', ...
                    qs0, qs, sl0, sl, (qs + sl) - (qs0 + sl0));
                
                obj.take_snapshot(problem);
            end            
        end
        
        function on_models_updated(obj, problem)  
            
            if obj.event_level >= 3
                lm0 = obj.log_mprior;
                sl0 = obj.sum_loglik;
                
                lm = problem.log_mprior;
                sl = problem.sum_loglik;
                
                fprintf(obj.output_fid, ...
                    '\tmodels updated: log_mprior %g --> %g, sum_loglik %g --> %g [ch = %g]\n', ...
                    lm0, lm, sl0, sl, (lm + sl) - (lm0 + sl0));
                
                obj.take_snapshot(problem);                
            end            
        end
        
        function on_completed(obj, problem) 
            
            if obj.event_level >= 1
                if problem.rtstatus.converged
                    fprintf(obj.output_fid, 'The procedure converged.\n');
                else
                    fprintf(obj.output_fid, 'The procedure did NOT converge.\n');
                end
            end
        end
        
        
        function take_snapshot(obj, problem)
            
            obj.log_mprior = problem.log_mprior;
            obj.log_qprior = problem.log_qprior;
            obj.qscore = problem.qscore;
            obj.sum_loglik = problem.sum_loglik;
            obj.total_objv = problem.total_objv;
        end
    end
    
end