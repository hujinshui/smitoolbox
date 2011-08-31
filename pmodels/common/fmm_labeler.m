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
        loglik_fun;     % the function handle for evaluating log-likelihood
                
        K;              % the number of models
        N;              % the number of observations
        opchar;         % the char indicating the method to use
                        % ('m' - 'mapest', 's' - 'sample', 'o'-'optima')
    end
    
    %% Construction
    
    methods
        
        function obj = fmm_labeler(fun, painfo, obsinfo, K, N, method)
            % Constructs an instance of FMM-based label inferrer
            %
            %   obj = fmm_labeler(f, painfo, obsinfo, K, N, method)
            %
            %   Input arguments:
            %       - fun:      the function handle to evaluate log
            %                   likelihood values pairwisely.
            %
            %       - painfo:   the struct of params information
            %
            %       - obsinfo:  the struct of observation information
            %
            %       - K:        the number of params to estimate
            %                   One can set K = 0 to indicate that the
            %                   such number is not fixed.
            %
            %       - N:        the number of observed samples.
            %                   One can set N = 0 to indicate that it is
            %                   not fixed. 
            %
            %       - method:   It can be either 'mapest' or 'sample',
            %
            %                   - 'mapest': the output is the MAP
            %                               estimation of the posterior
            %                               distribution of the labels
            %
            %                   - 'sample': the output is a sample from
            %                               posterior distribution of
            %                               labels.
            %
            %                   - 'optima': use the optimal label.
            %
            
            % verify input
            
            if ~isa(fun, 'function_handle')
                error('fmm_labeler:invalidarg', ...
                    'fun should be a function handle.');
            end
            
            if ~( isstruct(painfo) && ...
                    all(isfield(painfo, {'name', 'type', 'size'})) )
                error('fmm_labeler:invalidarg', ...
                    'painfo should be a variable-info struct.');
            end
            
            if ~( isstruct(obsinfo) && ...
                    all(isfield(obsinfo, {'name', 'type', 'size'})) )
                error('fmm_labeler:invalidarg', ...
                    'painfo should be a variable-info struct.');
            end
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
                error('fmm_labeler:invalidarg', ...
                    'K should be a non-negative integer scalar');
            end
            
            if ~(isnumeric(N) && isscalar(N) && N == fix(N) && N >= 0)
                error('fmm_labeler:invalidarg', ...
                    'N should be an integer scalar with N >= -1.');
            end
            
            if ~ischar(method)
                error('fmm_labeler:invalidarg', ...
                    'method should be a char string.');
            end
            
            switch lower(method)
                case 'mapest'
                    opc = 'm';
                case 'sample'
                    opc = 's';
                otherwise
                    error('fmm_labeler:invalidarg', ...
                        'The method is invalid.');
            end
            
            % create slots
            
            np = numel(painfo);
            m = numel(obsinfo);
            
            inlets = repmat(struct('name',[],'type',[],'size',[]), ...
                np+m+1, 1);
            
            for i = 1 : np
                inlets(i).name = painfo(i).name;
                inlets(i).type = painfo(i).type;
                inlets(i).size = [painfo(i).size, K];
            end
            
            for i = 1 : m
                inlets(np+i).name = obsinfo(i).name;
                inlets(np+i).type = obsinfo(i).type;
                inlets(np+i).size = [obsinfo(i).size, N];
            end
            
            inlets(np+m+1).name = 'LabelPrior';
            inlets(np+m+1).type = 'double';
            inlets(np+m+1).size = [K, 1];
            
            outlets(1).name = 'Labels';
            outlets(1).type = 'double';
            if opc == 'm'
                outlets(1).size = [K, N];
            else % == 's' or 'o'
                outlets(1).size = [1, N];
            end
            
            outlets(2).name = 'LogLik';
            outlets(2).type = 'double';
            outlets(2).size = [K, N];
            
            obj = obj@smi_func(inlets, outlets);
            
            % set fields
            
            obj.loglik_fun = fun;
            obj.K = K;
            obj.N = N;
            obj.opchar = opc;
                        
        end
        
        
        %% Implementation

        function tf = test_slots(obj, inflags, ~) %#ok<MANU>
            % Test whether a given input/output pattern is acceptable
                        
            tf = all(inflags(1:end-1));            
        end

                
        function [Z, LLik] = evaluate(obj, ~, varargin) 
            % Evaluate the function            
                        
            % get prior
            
            pri = varargin{end};
            
            % evaluate log-likelihood
            
            LLik = obj.loglik_fun(varargin{1:end-1});
            
            % evaluate Z
                       
            if isempty(pri)
                E = LLik;
            else
                E = bsxfun(@plus, log(pri), LLik);
            end
            
            opc = obj.opchar;
            
            if opc == 'o'
                [~, Z] = max(E, [], 1);
            else            
                Q = nrmexp(E, 1);
            
                if opc == 's'
                    Z = ddsample(Q, 1);
                else
                    Z = Q;
                end 
            end
        end
        
        
    end
    
    
end
