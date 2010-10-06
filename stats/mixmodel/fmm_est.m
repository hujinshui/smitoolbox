classdef fmm_est
    % The class for estimating a finite mixture model
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
        
    properties
        qlabeler;       % the labeler for inferring posterior label distribution
        model_est;      % the model estimator
        
        maxiter = 200;      % the maximum number of iterations
        tolerance = 1e-6;   % the tolerance of change of objective at convergence
        event_level = 3;    % the level of event notification
                            % 0 - no event is raised
                            % 1 - only procedure-level event is raised
                            % 2 - up to iteration-level event is raised
                            % 3 - up to step-level event is raised                            
    end
    
    methods
        function obj = set.qlabeler(obj, v)
            assert(isa(v, 'labeler'), 'fmm_est:invalidprop', ...
                'qlabeler should be an object of class labeler.');
            obj.qlabeler = v;
        end        
        
        function obj = set.model_est(obj, v)
            assert(isa(v, 'smodelest'), 'fmm_est:invalidprop', ...
                'model_est should be an object of class smodelest.');
            obj.model_est = v;
        end
        
        function obj = set.maxiter(obj, v)
            assert(isnumeric(v) && isscalar(v) && v >= 1, ...
                'fmm_est:invalidprop', ...
                'maxiter should be a positive number no less than 1.');
            obj.maxiter = v;
        end
        
        function obj = set.tolerance(obj, v)
            assert(isfloat(v) && isreal(v) && isscalar(v) && v > 0, ...
                'fmm_est:invalidprop', ...
                'maxiter should be a positive real scalar.');
            obj.tolerance = v;
        end
        
        function obj = set.event_level(obj, v)
            assert(isnumeric(v) && isscalar(v) && v >= 0, ...
                'fmm_est:invalidprop', ...
                'event_level should be a positive number.');
            obj.event_level = v;
        end
    end
    
    
    methods
        
        function Pb = accept(obj, X, initQ)
            % Check the validity of input samples and form an estimation
            % problem, and initialize the solution
            %
            %   Pb = obj.accept(X, initQ);
            %
            %   In the input, X is the observation in a form that can 
            %   be recognized by model_est, and initQ is the initial
            %   guess of the label assignment.
            %
            
            Pb = fmm_est_pb(obj, X, initQ);
        end        
        
    end
            
end
