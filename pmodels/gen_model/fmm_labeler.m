classdef fmm_labeler < smi_func
    % An smi function for inferring labels based on finite mixture model
    %
    %   Suppose there are 
    %       np parameters of each model
    %       m  products of each observation
    %
    %   (np + m + 1) input slots:
    %       1 : np          - model parameters
    %       (np+1) : (np+m) - observed products
    %       np+m+1          - prior of labels (optional)
    %
    %   2 output slots
    %       1 - labels (for sampling) or label distributions (otherwise)
    %       2 - pairwise log-likelihood of observations with all models
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')
        model;  % the generative model object (gen_model_base)
        
        K;              % the number of models
        N;              % the number of observations
        is_sampler;     % whether it is for sampling
    end
    
    %% Construction
    
    methods
        
        function obj = fmm_labeler(model, K, N, op)
            % Constructs an instance of FMM-based label inferrer
            %
            %   obj = fmm_labeler(model, K, N, method)
            %
            %   Input arguments:
            %       - model:    the generative model (gen_model_base)
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
            %                               estimation of the posterior
            %                               distribution of the labels
            %
            %                   - 'sample': the output is a sample from
            %                               posterior distribution of
            %                               labels.
            %
            
            % verify input
            
            if ~isa(model, 'gen_model_base')
                error('fmm_labeler:invalidarg', ...
                    'model should be an object of class gen_model_base.');
            end
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 0)
                error('fmm_labeler:invalidarg', ...
                    'K should be a non-negative integer scalar.');
            end
                        
            if ~(isnumeric(N) && isscalar(N) && N == fix(N) && N >= 0)
                error('fmm_labeler:invalidarg', ...
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
            
            np = numel(model.params_info);
            m = numel(model.products_info);
                                    
            inlets = repmat(struct('name', [], 'type', [], 'size', []), ...
                np + m + 1, 1);
            
            painfo = model.params_info;
            for i = 1 : np
                inlets(i).name = painfo(i).name;
                inlets(i).type = painfo(i).type;
                inlets(i).size = [painfo(i).size, K];
            end
            
            prinfo = model.products_info;
            for i = 1 : m
                inlets(np+i).name = prinfo(i).name;
                inlets(np+i).type = prinfo(i).type;
                inlets(np+i).size = [prinfo(i).size, N];
            end
            
            inlets(np+m+1).name = 'LabelPrior';
            inlets(np+m+1).type = 'double';
            inlets(np+m+1).size = K;
            
            if is_samp
                outlets(1).name = 'Labels';
                outlets(1).type = 'numeric';
                outlets(1).size = [1, N];
            else
                outlets(1).name = 'LabelDistr';
                outlets(1).type = 'double';
                outlets(1).size = [K, N];
            end
            
            outlets(2).name = 'LogLik';
            outlets(2).type = 'double';
            outlets(2).size = [K, N];
            
            % make object
            
            obj = obj@smi_func(inlets, outlets);
            
            obj.model = model;
            obj.K = K;
            obj.N = N;
            obj.is_sampler = is_samp;            
        end
        
        
        %% Implementation

        function tf = test_slots(obj, inflags, outflags) %#ok<INUSD>
            % Test whether a given input/output pattern is acceptable
            
            np = obj.model.num_params;
            m = obj.model.num_products;            
            tf = all(inflags(1:np+m));            
        end

                
        function [Z, LLik] = evaluate(obj, outflags, varargin) %#ok<INUSL>
            % Evaluate the function            
            
            % get inputs
            
            mdl = obj.model;
            
            np = mdl.num_params;
            m = mdl.num_products;
            
            params = varargin(1:np);
            obs = varargin(np+1:np+m);
            pri = varargin{np+m+1};
            
            % evaluate log-likelihood
            
            LLik = mdl.loglik(params, obs);
            
            % evaluate Z
                       
            if isempty(pri)
                E = LLik;
            else
                E = bsxfun(@plus, pri, LLik);
            end
            
            Q = nrmexp(E, 1);
            
            if obj.is_sampler
                Z = ddsample(Q, 1);
            else
                Z = Q;
            end                
        end
        
        
    end
    
    
end
