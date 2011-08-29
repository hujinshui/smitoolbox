classdef fmm_param_inferrer < smi_func
    % An smi function for inferring FMM model parameters
    %
    % Suppose there are nh distinct hyper-parameters.
    %   
    % (nh + 2) input slots:
    %   1 : nh - hyper parameters, 
    %   (nh+1) - observed data,
    %   (nh+2) - labels / label distribution
    % 
    % 1 output slot: to output inferred parameter    
    %
        
    % Created by Dahua Lin, on Aug 26, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')        
        model;      % the generative model (gen_model)        
        
        K;          % the number of parameters to estimate
        N;          % the number of observed samples               
        is_sampler;     % whether this object is an sampler
    end
    
    
    methods
        
        %% Construction
        
        function obj = fmm_param_inferrer(model, K, N, op)
            % Constructs a parameter estimator
            %
            %   obj = fmm_param_inferrer(model, K, op);
            %
            %   Input arguments:
            %       - model:    the generative model
            %
            %       - K:        the number of params to estimate
            %                   One can set K = 0 to indicate that the
            %                   such number is not fixed.
            %
            %       - N:        the number of observed samples.
            %                   One can set N = 0 to indicate that it is
            %                   not fixed. 
            %
            %       - op:       It can be either 'mapest' or 'sample',
            %
            %                   - 'mapest': the output is the MAP
            %                               estimation
            %                               (nh+2)-th input is posterior
            %                               probabilities of labels,
            %                               of size K x n.
            %
            %                   - 'sample': the output is a sample from
            %                               the posterior distribution.
            %                               (nh+2)-th input is a vector
            %                               of labels in [1, K].
            %
            
            % verify input
            
            if ~isa(model, 'gen_model')
                error('fmm_param_inferrer:invalidarg', ...
                    'model should be an object of class gen_model.');
            end
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 0)
                error('fmm_param_inferrer:invalidarg', ...
                    'K should be a non-negative integer scalar.');
            end
                        
            if ~(isnumeric(N) && isscalar(N) && N == fix(N) && N >= 0)
                error('fmm_param_inferrer:invalidarg', ...
                    'N should be a non-negative integer scalar.');
            end
            
            if strcmpi(op, 'mapest')
                is_samp = false;
            elseif strcmpi(op, 'sample')
                is_samp = true;
            else
                error('fmm_param_inferrer:invalidarg', ...
                    'Invalid op string.');
            end  
            
            % fill in slot information
            
            nh = model.num_hyper_params;
            
            inlets = repmat(struct('name', [], 'type', [], 'size', []), ...
                nh + 2, 1);
            
            hpinfo = model.hyper_params_info;
            for i = 1 : nh                
                inlets(i).name = hpinfo(i).name;
                inlets(i).type = hpinfo(i).type;
                inlets(i).size = hpinfo(i).size;
            end
            
            prinfo = model.products_info;
            inlets(nh+1).name = prinfo.name;
            inlets(nh+1).type = prinfo.type;
            inlets(nh+1).size = [prinfo.size, N];
            
            if is_samp
                inlets(nh+2).name = 'Labels';
                inlets(nh+2).type = 'numeric';
                inlets(nh+2).size = [1, N];                                
            else 
                inlets(nh+2).name = 'LabelDistr';
                inlets(nh+2).type = 'double';
                inlets(nh+2).size = [K, N];
            end
                           
            painfo = model.params_info;
            outlets.name = painfo.name;
            outlets.type = painfo.type;
            outlets.size = [painfo.size, K];            
            
            % make object
            
            obj = obj@smi_func(inlets, outlets);
            
            obj.model = model;
            obj.K = K;
            obj.N = N;            
            obj.is_sampler = is_samp;
        end
        
        
        %% Implementation

        function tf = test_slots(obj, inflags, outflags) %#ok<MANU,INUSD>
            % Test whether a given input/output pattern is acceptable
            %
            
            tf = all(inflags);            
        end

                
        function theta = evaluate(obj, outflags, varargin) %#ok<INUSL>
            % Evaluate the function
            %
            
            nh = obj.model.num_hyper_params;
            
            hparams = varargin(1:nh);
            obs = varargin{nh+1};
            Z = varargin{nh+2};
            
            if obj.is_samp;
                g = intgroup(obj.K, Z);
                theta = obj.model.sample_params(obs, g, hparams);
            else                
                theta = obj.model.mapest_params(obs, Z, hparams);
            end
        end
        
        
    end
    
    
end

