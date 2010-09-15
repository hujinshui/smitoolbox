classdef rfmm_est
    % The class for robust finite mixture model estimation
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
        
    properties
        qlabeler;       % the labeler for inferring posterior label distribution
        glabeler;       % the labeler for inferring inlier probabilities
        model_est;      % the model estimator
        lp_outlier;     % the log-probability of outliers
        
        maxiter = 200;      % the maximum number of iterations
        tolerance = 1e-6;   % the tolerance of change of objective at convergence
        event_level = 4;    % the level of event notification
                            % 0 - no event is raised
                            % 1 - only procedure-level event is raised
                            % 2 - only iteration-level event is raised
                            % 3 - only step-level event is raised                            
    end
    
    methods
        function obj = set.qlabeler(obj, v)
            assert(isa(v, 'labeler'), 'rfmm_est:invalidprop', ...
                'qlabeler should be an object of class labeler.');
            obj.qlabeler = v;
        end
        
        function obj = set.glabeler(obj, v)
            assert(isa(v, 'labeler'), 'rfmm_est:invalidprop', ...
                'glabeler should be an object of class labeler.');
            obj.glabeler = v;
        end
        
        function obj = set.model_est(obj, v)
            assert(isa(v, 'smodelest'), 'rfmm_est:invalidprop', ...
                'model_est should be an object of class smodelest.');
            obj.model_est = v;
        end         
        
        function obj = set.maxiter(obj, v)
            assert(isnumeric(v) && isscalar(v) && v >= 1, ...
                'rfmm_est:invalidprop', ...
                'maxiter should be a positive number no less than 1.');
            obj.maxiter = v;
        end
        
        function obj = set.tolerance(obj, v)
            assert(isfloat(v) && isreal(v) && isscalar(v) && v > 0, ...
                'rfmm_est:invalidprop', ...
                'maxiter should be a positive real scalar.');
            obj.tolerance = v;
        end
        
        function obj = set.event_level(obj, v)
            assert(isnumeric(v) && isscalar(v) && v >= 0, ...
                'rfmm_est:invalidprop', ...
                'event_level should be a positive number.');
            obj.event_level = v;
        end
    end
    
    
    methods
        
        function Pb = accept(obj, X, initQ, initG)
            % Check the validity of input samples and form an estimation
            % problem, and initialize the solution
            %
            %   Pb = obj.accept(X, initQ, initG);
            %
            %   In the input, X is the observation in a form that can 
            %   be recognized by model_est, and initQ is the initial
            %   guess of the label assignment, initG is the initial
            %   guess of the vector of inlier probabilities.  If all
            %   samples have the same confidence of being inlier, 
            %   then initG can be input as a scalar.
            %
            
            Pb = rfmm_est_pb(obj, X, initQ, initG);
        end        
        
    end
            
end
