classdef fmm
    % The class to represent finite mixture model
    %
    
    % Created by Dahua Lin, on Nov 14, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        gm;         % the generative model
        params;     % the model parameters (pd x K)
        
        ncomps;     % the number of components (K)
        priws;      % the prior weights of different components (1 x K)
        
    end
    
    
    methods
        
        function obj = fmm(gm, params, ws)
            % Constructs a finite mixture model
            %
            %   obj = fmm(gm, params, ws);
            %       constructs a finite mixture model.
            %
            %       Inputs:
            %       - gm:       the underlying generative model object
            %                   (such as gaussgm)
            %       - params:   the parameters of all components (pd x K)
            %       - ws:       the prior weights of all components (1 x K)
            %                   Note that sum(ws) should equal 1.
            %   
            
            if ~isobject(gm)
                error('fmm:invalidarg', 'gm should be an object.');
            end
            
            K = gm.check_params(params);
            if K <= 0
                error('fmm:invalidarg', 'the params are invalid.');
            end
            
            if ~(isfloat(ws) && ndims(ws) == 2 && size(ws,1) == 1)
                error('fmm:invalidarg', 'ws should be a row vector.');
            end
                       
            if K ~= size(ws, 2)
                error('fmm:invalidarg', ...
                    'The length of ws should equal the number of components.');
            end
            
            obj.gm = gm;
            obj.params = params;
            obj.ncomps = K;
            obj.priws = ws;
        end
        
        
        function lp = logpdf(obj, X)
            % Compute the log pdf (of the mixed distribution) at given samples
            %
            %   lp = obj.logpdf(X);
            %       
            %       Suppose there are n samples, then lp will be 
            %       a row vector of size 1 x n.
            %
                                 
            LL = obj.gm.loglik(obj.params, X);            
            maxL = max(LL, [], 1);
            DL = bsxfun(@minus, LL, maxL);
            
            lp = maxL + log(obj.priws * exp(DL));                        
        end
    
    end
    
end