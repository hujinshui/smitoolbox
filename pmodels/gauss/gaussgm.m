classdef gaussgm < genmodel_base
    % The class to implement a basic Gaussian generative model
    %
    %   x ~ N(mu, Sigma);
    %
    %   Two parameters:
    %   - mu:       the mean vector
    %   - Sigma:    the covariance 
    %
    %   They together are encapsulated in a gaussd struct.
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        dim;                % the space dimension   
        cov_form;           % the form of covariance
        tied_cov = false;   % whether the covariance is tied
    end
    
    %% constructor
    
    methods
        
        function model = gaussgm(d, cf, op)
            % Constructs a Gaussian generative model
            %
            %   model = gaussgm(d, cf);
            %       creates a Gaussian generative model of specified
            %       dimension and covariance form.
            %
            %       Inputs:
            %       - d:        the space dimension
            %       - cf:       the form of covariance:
            %                   's':    isotropic covariance
            %                   'd':    diagonal covariance
            %                   'f':    full-form covariance
            %
            %   model = gaussgm(d, cf, 'tied-cov');
            %       creates a Gaussian generative model where the
            %       covariance of all component models are tied.
            %       
           
            if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1)
                error('gaussgm:invalidarg', 'd should be a positive integer.');
            end
            
            if ~(ischar(cf) && isscalar(cf) && any(cf == 'sdf'))
                error('gaussgm:invalidarg', ...
                    'cf should be either of ''s'', ''d'', or ''f''.');
            end
            
            model.dim = d;
            model.cov_form = cf;
            
            if nargin >= 3
                if ~strcmpi(op, 'tied-cov')
                    error('gaussgm:invalidarg', ...
                        'The third argument is invalid.');
                end
                model.tied_cov = true;
            end            
        end
    end
    
    
    %% Query and Evaluation
    
    methods
        
        function n = query_obs(model, X)
            % Get the number of samples in the input 
            %
            %   n = model.query_obs(X);
            %
            
            d = model.dim;
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == d)
                error('gaussgm:invalidarg', ...
                    'The sample matrix should be a real matrix with d rows.');
            end            
            n = size(X, 2);            
        end
        
        function n = query_params(model, G)
            % Get the number of parameters in the input
            %
            %   n = model.query_params(G);
            %
            
            d = model.dim;
            if ~(is_gaussd(G) && G.ty == 'm' && G.d == d)
                error('gaussgm:invalidarg', ...
                    'The parameters G should be a gaussd struct with G.d == d.');
            end 
            n = G.n;
        end
                
        function LL = loglik(model, G, X)
            % Evaluate the log-likelihood values at given samples
            %
            %   LL = loglik(model, G, X);
            %
            %       evaluates the log-likelihood at the samples given
            %       in X, with respect to the Gaussian distributions
            %       represented by G.
            %
            
            d = model.dim;
            if ~(is_gaussd(G) && G.d == d)                                
                error('gaussgm:invalidarg', ...
                    'The parameters G should be a gaussd struct with G.d == d.');
            end            
            LL = gaussd_logpdf(G, X);
        end
        
    end
    
    
    %% Estimation
    
    methods
        
        function G = mle(model, X, Z)
            % Performs maximum likelihood estimation of the parameters
            %
            %   G = model.mle(X, Z);
            %
            %       performs maximum likelihood estimation based on
            %       given (weighted) set of data
            %
            
            n = model.query_obs(X);
            [zty, K] = verify_Zarg(Z, n);
            
            if zty <= 1
                w = Z;
            else
                w = zeros(K, n);
                for k = 1 : K
                    w(k, Z{k}) = 1;
                end
            end
            
            cf = model.cov_form;
            tie_c = model.tied_cov;
            
            G = gaussd_mle(X, w, cf, tie_c);            
        end
        
        
        function capture(model, X, Z) %#ok<INUSD,MANU>
            error('gaussgm:notsupported', ...
                'The capture method is not supported by gaussgm');
        end
                        
    end        
    
end


