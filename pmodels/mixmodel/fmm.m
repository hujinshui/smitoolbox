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
            
            if ~(isfloat(params) && ndims(params) == 2)
                error('fmm:invalidarg', 'params should be a numeric matrix.');
            end
            
            if ~(isfloat(ws) && ndims(ws) == 2 && size(ws,1) == 1)
                error('fmm:invalidarg', 'ws should be a row vector.');
            end
            
            K = size(params, 2);
            if K ~= size(ws, 2)
                error('fmm:invalidarg', ...
                    'The size of params and that of ws are inconsistent.');
            end
            
            obj.gm = gm;
            obj.params = params;
            obj.ncomps = K;
            obj.priws = ws;            
        end
    
    end
    
end