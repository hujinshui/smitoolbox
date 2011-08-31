classdef gaussgm_mean_inferrer < smi_func
    % The SMI function for inferring mean parameter for Gaussian model
    %
    
    % Created by Dahua Lin, on Aug 30, 2011
    %
    
    %% Properties
    
    properties
        model;      % the generative model (of class gaussgm)
        
        K;          % the number of parameters to infer
        N;          % the number of observed samples
        opchar;     % the char indicating the method of doing inference
                    % 'm':  performs MAP (or ML) estimation
                    % 's':  samples from posterior distribution
    end
    
    %% Construction
    
    methods
        
        function obj = gaussgm_mean_inferrer(model, K, N, method)
            % Constructs a gaussgm_mean_inferrer object
            %
            %   obj = gaussgm_mean_inferrer(model, K, N, method);
            %
            %   Input arguments:
            %   - model:    the underlying generative model, which 
            %               should be of class gaussgm
            %
            %   - K:        the number of parameters (mean vectors) to
            %               infer. Set K = -1 when this number is not
            %               fixed.
            %
            %   - N:        the number of observed samples. Set N = -1
            %               when this number is not fixed.
            %
            %   - method:   the inference method to use
            %               'mapest':   performs MAP (or ML) estimation
            %               'sample':   samples from the posterior
            %
            
            % verify inputs
            
            if ~isa(model, 'gaussgm')
                error('gaussgm_mean_inferrer:invalidarg', ...
                    'model should be an object of gaussgm.');
            end
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= -1)
                error('gaussgm_mean_inferrer:invalidarg', ...
                    'K should be an integer scalar with K >= -1.');
            end
            
            if ~(isnumeric(N) && isscalar(N) && N == fix(N) && N >= 0)
                error('gaussgm_mean_inferrer:invalidarg', ...
                    'N should be an integer scalar with N >= -1.');
            end
            
            if ~ischar(method)
                error('gaussgm_mean_inferrer:invalidarg', ...
                    'method should be a char string.');
            end
            
            switch lower(method)
                case 'mapest'
                    opc = 'm';
                case 'sample'
                    opc = 's';
                otherwise
                    error('gaussgm_mean_inferrer:invalidarg', ...
                        'The method is invalid.');
            end
            
            % create object
                                    
            inlets = struct( ...
                'name', {'x', 'mu', 'z', 'Cx'}', ...
                'type', [], 'size', []);
            
            inlets(1).type = 'double';
            inlets(1).size = [model.xdim, N];
            
            inlets(2).type = 'double';
            inlets(2).size = model.dim0;
            
            if opc == 'm'
                inlets(3).type = 'double';
                inlets(3).size = [K, N];
            else % opc == 's'
                inlets(3).type = 'numeric';
                inlets(3).size = [1, N];
            end
            
            inlets(4).type = 'struct';
            inlets(4).size = 1;
            
            outlet.name = 'u';
            outlet.type = 'double';
            outlet.size = [model.udim, K];
            
            obj = obj@smi_func(inlets, outlet); 
            
            obj.model = model;
            obj.K = K;
            obj.N = N;
            obj.opchar = opc;
        end
        
    end
    
    
    %% Evaluation
    
    methods
        
        function tf = test_slots(obj, inflags, ~) 
            % Test whether the slots are properly connected
            
            tf_cx = obj.model.fixed_Cx || inflags(4);
            
            if obj.K == 1
                tf = inflags(1) && tf_cx;
            else
                tf = inflags(1) && inflags(3) && tf_cx;
            end
        end
        
        
        function U = evaluate(obj, ~, X, Mu, Z, Cx)
            % Performs the inference
            %
            %   U = obj.evaluate(outflags, X, Mu, Z);
            %
            %   Inputs:
            %   - X:        observed samples [xd x n]
            %   - Mu:       prior means [d0 x K] (optional)
            %   - Z:        sample labels [1 x n or K x n] (optional)
            %
            %   Outputs:
            %   - U:        the inferred/sampled mean parameter
            %
            
            opc = obj.opchar;
            if opc == 'm'
                U = obj.model.mapest_u(X, Cx, Mu, Z);
            elseif opc == 's'
                g = intgroup(obj.K, Z);
                U = obj.model.sample_u(X, Cx, Mu, g);
            end
        end
    
    end
        
end


