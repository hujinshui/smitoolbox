classdef ppca_gm < genmodel_base
    % The class implementing PPCA in terms of genmodel_base
    %
    
    % Created by Dahua Lin, on Nov 6, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')        
        dim;        % the observed space dimension
        ldim;       % the latent space dimension                
    end
    
    
    methods
        
        function model = ppca_gm(d, q)
            % Construct a PPCA generative model
            %
            %   model = ppca_gm(d, q);
            %       constructs a PPCA generative model.
            %       
            %       - d:    the observed space dimension
            %       - q:    the latent space dimension
            %
            
            model.dim = d;
            model.ldim = q;            
        end
        
        
        function tf = is_supported_optype(model, optype) %#ok<MANU>
            % Tests whether a particular operation type is supported
            %
            %   returns true only when optype is 'varinfer'
            
            tf = ischar(optype) && strcmpi(optype, 'varinfer');
        end
        
        
        function n = get_num_observations(model, obs) 
            % Gets the number of samples in the observation set
            %
            %   n = model.get_num_observations(obs);
            %       
            
            if ~(isfloat(obs) && ndims(obs) == 2 && isreal(obs) && ~issparse(obs))
                error('genmodel_base:invalidarg', ...
                    'obs should be a non-sparse real matrix.');
            end            
            
            if size(obs,1) ~= model.dim
                error('genmodel_base:invalidarg', ...
                    'The dimension of obs is invalid.');
            end
            
            n = size(obs, 2);
        end
        
        
        function tf = is_valid_prior(model, pri) %#ok<MANU>
            % Tests whether a given prior is a valid prior model
            %
            %   tf = model.is_valid_prior(pri)
            %
            %   Currently, it only accepts empty prior.
            %
            
            tf = isempty(pri);
        end
        
        
        function [params, aux] = posterior_params(model, ~, X, aux, Z, optype)
            % Estimates/samples the parameters conditoned on given observations        
            %
            %   [params, aux] = model.posterior_params([], X, aux, Z, 'varinfer');
            %       estimates/infer K parameters given grouped/weighted data
            %       set.
            %
            %       Input Arguments:
            %       - obs:      the observation data set
            %       - Z:        A K x n matrix of weights.
            %
            %       Output Arguments:
            %       - params:   a cell array, each cell is a probpca 
            %                   object.                                          
            %
            %       - aux:      an updated auxiliary structure, 
            %                   which is [].        
            %
            
            if ~strcmpi(optype, 'varinfer')
                error('ppca_gm:invalidarg', 'Unsupported optype');
            end
            
            q = model.ldim;
            K = size(Z, 1);
            
            params = cell(1, K);
                        
            for k = 1 : K
                z = Z(k, :);
                params{k} = probpca.mle(X, z, q);
            end
        
        end
                
        
        function Lliks = evaluate_logliks(model, params, obs, ~, I) %#ok<MANU>
            % Evaluates the logarithm of likelihood of samples
            %
            %   Lliks = model.evaluate_logliks(params, obs, []);
            %   Lliks = model.evaluate_logliks(params, obs, [], I);
            %
            %       evaluates the log-likelihood values of the observations
            %       in obs with respect to the given set of parameters or
            %       the given atom.
            %
            %       If params has K parameters, Lliks is a K x n matrix.
            %
            
            if nargin < 5
                X = obs;
            else
                X = obs(:, I);
            end
            
            K = numel(params);
            n = size(X, 2);
            Lliks = zeros(K, n);
            
            for k = 1 : K
                pm = params{k};
                Lliks(k,:) = logpdf(pm, X);
            end            
        end
        
        
        
        function lpri = evaluate_logpri(model, ~, ~, ~) %#ok<MANU>
            % Evaluate the total log-prior of a given set of parameters
            %
            %   returns 0
        
            lpri = 0;
        end
        
    end
    
    
end
