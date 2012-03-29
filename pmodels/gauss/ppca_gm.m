classdef ppca_gm < genmodel_base
    % The class implementing PPCA as a genmodel_base subclass
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
        
        
        function n = query_obs(model, obs) 
            % Gets the number of samples 
            %
            %   n = model.query_obs(obs);
            %       
            
            if ~(isfloat(obs) && ndims(obs) == 2 && isreal(obs))
                error('ppca_gm:invalidarg', ...
                    'obs should be a non-sparse real matrix.');
            end            
            
            if size(obs,1) ~= model.dim
                error('ppca_gm:invalidarg', ...
                    'The dimension of obs is invalid.');
            end
            
            n = size(obs, 2);
        end
        
        
        function n = query_params(model, Ms)
            % Gets the number of models
            %
            %   n = model.query_params(Ms);
            %       
            
            d = model.dim;
            q = model.ldim;
            
            if ~(isstruct(Ms) && is_ppca(Ms(1)) && Ms(1).d == d && Ms(1).q == q)
                error('ppca_gm:invalidarg', ...
                    'The model (params) are invalid.');
            end
            
            n = numel(Ms);
        end
                        
        
        function L = loglik(model, Ms, X) %#ok<MANU>
            % Evaluate the log-likelihood of all samples w.r.t all models
            %
            %   L = model.loglik(Ms, X);
            %
            
            K = numel(Ms);
            if K == 1
                L = ppca_logpdf(Ms, X);
            else
                L = zeros(K, size(X,2));
                for k = 1 : K
                    L(k,:) = ppca_logpdf(Ms(k), X);
                end
            end
        end
        
        
        function Ms = mle(model, X, W, I)
            % Performs maximum likelihood estimation of models
            %
            %   Ms = model.mle(X, W);
            %   Ms = model.mle(X, W, I);
            %
            
            if nargin < 3
                W = [];
            end
                        
            n = model.query_obs(X);
            
            if ~isempty(W)
                if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,1) == n)
                    error('ppca_gm:invalidarg', ...
                        'W should be a real matrix with n columns.');
                end
                K = size(W, 2);
            else
                W = [];
                K = 1;
            end
            
            if nargin >= 4
                X = X(:, I);
                if ~isempty(W)
                    W = W(I, :);
                end
            end
            
            q = model.ldim;
            if K == 1
                Ms = ppca_mle(X, W, q);
            else
                Ms = cell(1, K);
                for k = 1 : K
                    Ms{k} = ppca_mle(X, W(:,k), q);
                end
                Ms = vertcat(Ms{:});
            end
        end
        
        
        function capture(model, X, W) %#ok<INUSD,MANU>
            error('ppca_gm:notsupported', ...
                'The capture method is not supported by PPCA.');            
        end
        
    end    
    
end


