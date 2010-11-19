classdef gaussgm_gp
    % The class to represent a Gaussian generative model with Gauss prior
    %
    %   The model is formalized as follows:
    %
    %       theta ~ N(mu_pri, sigma_pri);  (prior of parameters)  
    %
    %       x ~ N(theta, sigma);  (generation of observation)
    %       or
    %       x ~ N(A * theta, sigma);
    %   
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Nov 13, 2010
    %       - Modified by Dahua Lin, on Nov 15, 2010
    %           - supports multi-prior.
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        pdim;       % the parameter dimension
        xdim;       % the observation dimension
                
        num_priors; % the number of prior distributions
        prior;      % the prior Gaussian model         
        A;          % the transform matrix 
                    % (can be empty, if the transform is identity)
                    
        sigma_x;    % the covariance of measurement noise   
        isigma_x;   % the inverse of sigma_x
        gnoise;     % the Gaussian noise distribution
        
        gmargin;    % the marginal Gaussian distribution
    end
    
    
    methods
        
        %% Construction and Basic manipulation
        
        function obj = gaussgm_gp(prior, A, sigma_x)
            % Constructs the Gaussian generative model
            %
            %   obj = gaussgm(prior, A, sigma);
            %
            %       constructs a Gaussian generative model as 
            %       formalized above.
            %       
            %       Inputs:
            %       - prior:    a gaussd object (prior.has_ip == true)
            %                   to represent the prior of parameters.
            %
            %       - A:        the generative transform, which 
            %                   can be either empty or a matrix
            %
            %       - sigma:    the covariance of measurement noise 
            %                   (can be a scalar, matrix, or matrix object).
            %
            
            if ~(isa(prior, 'gaussd') && prior.has_ip)
                error('gaussgm:invalidarg', ...
                    'prior should be a gaussd object with information parameters.');
            end
            
            pd = prior.dim;
            
            if isempty(A)
                A = [];
                xd = pd;
            else
                if ~(isfloat(A) && ndims(A) == 2 && size(A,2) == pd)
                    error('gaussgm:invalidarg', ...
                        'A should be a numeric matrix with size(A,2) == pdim.');
                end
                xd = size(A,1);
            end
            
            if isfloat(sigma_x)
                if isscalar(sigma_x)
                    sigma_x = udmat(xd, sigma_x);
                else
                    if ~isequal(size(sigma_x), [xd xd])
                        error('gaussgm:invalidarg', ...
                            'The size of sigma_x is invalid.');
                    end
                    sigma_x = gsymat(sigma_x);
                end
            elseif isobject(sigma_x)
                if ~(sigma_x.d == xd && sigma_x.n == 1)
                    error('gaussgm:invalidarg', 'sigma_x is invalid.');
                end
            else
                error('gaussgm:invalidarg', 'sigma_x is invalid.');
            end                        
                
            obj.pdim = pd;
            obj.xdim = xd;
            
            obj.prior = prior;
            obj.num_priors = prior.num;
            obj.A = A;
            obj.sigma_x = sigma_x;    
            obj.isigma_x = inv(sigma_x);
            obj.gnoise = gaussd.from_ip(0, obj.isigma_x, 0, 'mp');
            
            % compute marginal
            
            mu_pri = prior.mu;
            if isequal(mu_pri, 0)
                Amu = 0;
            else
                if isempty(A);
                    Amu = mu_pri;
                else
                    Amu = A * mu_pri;
                end
            end
            
            if isempty(A)
                mg_sigma = prior.sigma + sigma_x;
            else                                
                ASA = gsymat(qtrans(prior.sigma, A));                                    
                mg_sigma = ASA + sigma_x;
            end
            
            obj.gmargin = gaussd.from_mp(Amu, mg_sigma, 'ip');            
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
            if ~(isfloat(thetas) && ...
                    (isequal(thetas, 0) || ndims(thetas) == 2 && size(thetas,1) == pd))
                error('gaussgm:loglik:invalidarg', ...
                    'thetas should be a pdim x m numeric matrix.');
            end
            
            J = obj.isigma_x;
            
            A_ = obj.A;
            if isequal(thetas, 0)
                Jmu = 0;
            else
                if isempty(A_)
                    mu = thetas;
                else
                    mu = A_ * thetas;
                end
                Jmu = J * mu;
            end
                        
            if nargin < 3
                Gs = gaussd.from_ip(Jmu, J);
            else
                Gs = gaussd.from_ip(Jmu, J, [], ump);
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
                LP = pri.logpdf_map(thetas, pri_map);
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
        
        
        
        function P = get_posterior(obj, X, w, k, ump)
            % Get the posterior distribution of parameters
            %
            %   P = obj.get_posterior(X, w);
            %   P = obj.get_posterior(X, w, k);
            %   P = obj.get_posterior(X, w, k, 'mp');
            %       
            %       computes the posterior distribution of parameters
            %       based on observed data given in X, with respect to
            %       the k-th prior. (When multi-prior exist, k must be
            %       explicitly given).
            %
            %       The observed samples can be weighted by w. If there 
            %       are n samples, then X should be a matrix of size 
            %       xdim x n. 
            %
            %       In weighted case, w should be a row vector of size
            %       1 x n.            
            %
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == obj.xdim)
                error('gaussgm:invalidarg', ...
                    'X should be an d x n numeric matrix.');
            end
            
            if nargin < 3
                w = 1;
            end                      
            [ha, Ja] = gcondupdate(obj.A, X, obj.isigma_x, w);
            
            pri = obj.prior;
            if nargin < 4
                k = [];
            end                        
            if isempty(k)    
                if pri.num ~= 1
                    error('gaussgm:get_posterior:invalidarg', ...
                        'When there are multi priors, k must be specified.');
                end                
            end
            
            if nargin < 5
                P = pri.posterior(ha, Ja, k);                
            else
                P = pri.posterior(ha, Ja, k, ump);
            end            
        end
                
        
        function theta = estimate_map(obj, X, w, pri_map)
            % Performs MAP estimation of parameter
            %
            %   theta = obj.infer(X);
            %   theta = obj.infer(X, w);
            %   theta = obj.infer(X, w, pri_map);
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
            if nargin < 4
                pri_map = [];
            end
            
            K = size(w, 1);
            if K == 1            
                [ha, Ja] = gcondupdate(obj.A, X, obj.isigma_x, w);
                theta = obj.prior.pos_mean(ha, Ja, pri_map);
            else
                theta = zeros(obj.pdim, K);
                for k = 1 : K
                    [ha, Ja] = gcondupdate(obj.A, X, obj.isigma_x, w(k,:));
                    theta(:,k) = obj.prior.pos_mean(ha, Ja, pri_map);
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
            
            A_ = obj.A;
            if isempty(A_)
                mu = theta;
            else
                mu = A_ * theta;
            end

            e = obj.gnoise.sample(n);
            X = bsxfun(@plus, mu, e);            
            
        end
        
        
        function Y = pri_sample(obj, n, k)
            % Draw samples from the parameter priors
            %
            %   Y = obj.pri_sample(n, k);
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
                [ha, Ja] = gcondupdate(obj.A, X, obj.isigma_x, w);
                Y = obj.prior.pos_sample(ha, Ja, n, k);
            else
                Y = obj.prior.sample(n, k);
            end
        end
        
    end
end 

