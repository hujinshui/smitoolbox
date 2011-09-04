classdef gaussgm_gp
    % The class to represent a Gaussian generative model with Gauss prior
    %
    %   The model is formalized as follows:
    %
    %       theta ~ N(mu_pri, C_pri);  (prior of parameters)  
    %
    %       x ~ N(theta, C_x);  (generation of observation)
    %   
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Nov 13, 2010
    %       - Modified by Dahua Lin, on Nov 15, 2010
    %           - supports multi-prior.
    %       - Modified by Dahua Lin, on Aug 17, 2011
    %           - based on new gaussd object inferface.
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        pdim;       % the parameter dimension
        xdim;       % the observation dimension
                
        npri;       % the number of prior distributions
        prior;      % the prior Gaussian model         
                    
        gnoise;     % the Gaussian noise distribution        
        gmargin;    % the marginal Gaussian distribution
    end
    
    
    methods
        
        %% Construction and Basic manipulation
        
        function obj = gaussgm_gp(prior, Cx)
            % Constructs the Gaussian generative model
            %
            %   obj = gaussgm(prior, Cx);
            %
            %       constructs a Gaussian generative model as 
            %       formalized above.
            %       
            %       Inputs:
            %       - prior:    a gaussd object (prior.has_ip == true)
            %                   to represent the prior of parameters.           
            %
            %       - Cx:       the observation noise covariance  
            %                   (can be a scalar, a vector of diagonal
            %                    entries, or a full covariance matrix).
            %
            
            if ~(isa(prior, 'gaussd') && prior.has_ip)
                error('gaussgm_gp:invalidarg', ...
                    'prior should be a gaussd object with information parameters.');
            end
            
            pd = prior.dim;
            xd = pd;           
           
            if ~(isfloat(Cx) && ndims(Cx) == 2 && isreal(Cx))
                error('gaussgm_gp:invalidarg', ...
                    'Cx must be a real scalar/vector/matrix.');
            end
            
            if isscalar(Cx)
                cfx = 's';
            elseif size(Cx,1) == xd && size(Cx,2) == 1
                cfx = 'd';
            elseif size(Cx,1) == xd && size(Cx,2) == xd
                cfx = 'f';
            else
                error('gaussgm_gp:invalidarg', ...
                    'The size of Cx is invalid.');
            end
            
            % construct object
                                        
            obj.pdim = pd;
            obj.xdim = xd;
            
            obj.npri = prior.num;
            obj.prior = prior;

            obj.gnoise = gaussd.from_mp(cfx, zeros(xd,1), Cx, 'ip');
            
            % compute marginal distribution of x
                        
            [Cmg, mgf] = gmat_plus(prior.C, prior.cform, Cx, cfx);            
            obj.gmargin = gaussd.from_mp(mgf, prior.mu, Cmg, 'ip');            
        end
        
        
        function n = check_parameters(obj, params)
            % Check validity of parameters and return the number
            %
            %   n = obj.check_parameters(params);
            %       if params is a valid parameter matrix, it returns
            %       the number of params, otherwise it returns -1.
            %
            
            pd = obj.pdim;
            if isfloat(params) && ndims(params) == 2 && size(params,1) == pd
                n = size(params, 2);
            else
                n = -1;
            end
            
        end
        
        
        function n = check_observations(obj, X)
            % Check validity of observations and return the number
            %
            %   n = obj.check_observations(X);
            %       if X is a valid observation matrix, it returns
            %       the number of samples in X, otherwise it returns -1.
            %
            
            xd = obj.xdim;
            if isfloat(X) && ndims(X) == 2 && size(X,1) == xd
                n = size(X,2);
            else
                n = -1;
            end
            
        end

        
        function S = select_observations(obj, X, s) %#ok<MANU>
            % Select a subset of observations
            %
            %   S = obj.select_observations(X, s);
            %
            
            S = X(:, s);
        end
        
        
        function Gs = param_models(obj, thetas, ump)
            % Get Gaussian models with given parameters
            %
            %   Gs = obj.param_models(thetas);
            %       gets the Gaussian distribution models with given
            %       parameters
            %
            %   Gs = obj.param_models(thetas, 'mp');
            %       gets the Gaussian distribution models with given
            %       parameters (together with mean parameters)
            %
            
            pd = obj.pdim;
            if ~(isfloat(thetas) && ndims(thetas) == 2 && size(thetas,1) == pd)
                error('gaussgm_gp:param_models:invalidarg', ...
                    'thetas should be a pdim x m numeric matrix.');
            end
            
            cf = obj.gnoise.cform;
            J = obj.gnoise.J;
            h = J * thetas;
                                    
            if nargin < 3
                Gs = gaussd.from_ip(cf, h, J);
            else
                Gs = gaussd.from_ip(cf, h, J, [], ump);
            end
                
        end  
    end
        
                    
    methods
        
        %% Evaluation, Inference, Estimation & Sampling
                        
        function LP = logpri(obj, thetas, pri_map)
            % Compute the log-prior of parameters
            %
            %   LP = obj.logpri(thetas);
            %   LP = obj.logpri(thetas, pri_map);
            %
            
            pri = obj.prior;
            if nargin < 3
                if pri.num ~= 1
                    error('gaussgm_gp:invalidarg', ...
                        'pri_map is required when there are multi-priors.');
                end                
                LP = pri.logpdf(thetas);                 
            else
                LP = pri.logpdf_map(thetas, 1:size(thetas, 2), pri_map);
            end
            
        end            
        
        
        function LL = loglik(obj, thetas, X)
            % Compute log-likelihood
            %
            %   LL = obj.loglik(params, X);
            %       computes the log-likelihood at observed samples
            %       given by X with respect to the model parameters
            %       given by parameters.
            %
            %       Suppose there are m parameters, and n samples, then
            %       LL will be a matrix of size m x n.
            %
            
            Gs = obj.param_models(thetas);
            LL = Gs.logpdf(X);
        end
        
        
        function LM = logmargin(obj, X, k)
            % Compute log-marginal-likelihood
            %
            %   LM = obj.logmargin(X);
            %       computes the log marginal likelihood of the samples
            %       given in X with respect to all priors.
            %
            %   LM = obj.logmargin(X, k);
            %       computes the log marginal likelihood of the samples
            %       given in X with respect to the prior specified by k.
            %
            
            if nargin < 3
                LM = obj.gmargin.logpdf(X);
            else
                LM = obj.gmargin.logpdf(X, k);
            end
        end
                                
        
        function theta = estimate_map(obj, X, w)
            % Performs MAP estimation of parameter
            %
            %   theta = obj.infer(X);
            %   theta = obj.infer(X, w);
            %
            %       solves the maximum-a-posteriori estimation of
            %       parameters based on observed data given in X.            
            %       The observed samples can be weighted by w.
            %
            %       w can be a 1 x n row vector or a K x n matrix.                               
            %       In output, theta is a column vector of size pdim x 1,
            %       or size pdim x K.
            %
            %
                 
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == obj.xdim)
                error('gaussgm:invalidarg', ...
                    'X should be an d x n numeric matrix.');
            end
            
            if nargin < 3
                w = 1;
            end
            
            K = size(w, 1);
            Jx = obj.gnoise.J;
            cf = obj.gnoise.cform;
            
            if K == 1            
                [ha, Ja] = gcondupdate([], X, Jx, cf, w);
                theta = obj.prior.pos_mean(ha, Ja, cf);
            else
                theta = zeros(obj.pdim, K);
                for k = 1 : K
                    [ha, Ja] = gcondupdate([], X, Jx, cf, w(k,:));
                    theta(:,k) = obj.prior.pos_mean(ha, Ja, cf);
                end
            end            
        end                        
                       
        
        function X = sample(obj, theta, n)
            % Sample observations with given parameter
            %
            %   X = obj.sample(theta, n);
            %
            %       draw n samples from the likelihood model with
            %       parameter given by theta. 
            %
            
            pd = obj.pdim;
            if ~(isfloat(theta) && isequal(size(theta), [pd 1]))
                error('gaussgm_gp:sample:invalidarg', ...
                    'theta should be a column vector of size pd x 1.');
            end

            e = obj.gnoise.sample(n);
            X = bsxfun(@plus, theta, e);            
            
        end
        
        
        function Y = pri_sample(obj, n, k)
            % Draw samples from the parameter priors
            %
            %   Y = obj.pri_sample(n);
            %
            
            if nargin < 3
                k = [];
            end
            
            Y = obj.prior.sample(n, k);
        end
        
    
        function Y = pos_sample(obj, X, w, n, k)
            % Sampling from posterior distribution of parameters
            %
            %   Y = obj.pos_sample(X, w, n);
            %   Y = obj.pos_sample(X, w, n, k);
            %
            %       draw n samples from the posterior distribution
            %       of parameters given (weighted) observations in X.
            %       One can input w as a scalar, when all observations
            %       have the same weight.
            %
            %
            
            % verify arguments
            if ~(isempty(X) || (isfloat(X) && ndims(X) == 2 && size(X,1) == obj.xdim))
                error('gaussgm:invalidarg', ...
                    'X should be empty or an d x n numeric matrix.');
            end
            if nargin < 3; w = 1; end
            if nargin < 4; n = 1; end
            if nargin < 5; k = []; end            
            
            % do sampling
            
            if ~isempty(X)
                Jx = obj.gnoise.J;
                cf = obj.gnoise.cform;
                [ha, Ja] = gcondupdate([], X, Jx, cf, w);
                Y = obj.prior.pos_sample(ha, Ja, cf, n, k);
            else
                Y = obj.prior.sample(n, k);
            end
        end
        
    end
end 

