classdef gaussgm_gp
    % The class to represent a Gaussian generative model with Gauss prior
    %
    %   The model is formalized as follows:
    %
    %       theta ~ N(mu0, sigma0);  (prior of parameters)  
    %
    %       x ~ N(theta, sigma);  (generation of observation)
    %       or
    %       x ~ N(A * theta, sigma);
    %
    
    % Created by Dahua Lin, on Nov 13, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        pdim;       % the parameter dimension
        xdim;       % the observation dimension
        
        prior;      % the prior Gaussian model
        
        A;          % the transform matrix 
                    % (can be empty, if the transform is identity)
                    
        isigma;     % the inverse covariance of measurement noise    
        gnoise;     % the Gaussian noise distribution
    end
    
    
    methods
        
        function obj = gaussgm_gp(prior, A, isigma)
            % Constructs the Gaussian generative model
            %
            %   obj = gaussgm(prior, A, isigma);
            %       constructs a Gaussian generative model as 
            %       formalized above.
            %       
            %       Inputs:
            %       - prior:    a gaussd object (prior.has_ip == true)
            %                   to represent the prior of parameters
            %       - A:        the generative transform, which 
            %                   can be either empty or a matrix
            %       - isigma:   the inverse covariance of measurement
            %                   noise (can be a scalar, matrix, or
            %                   matrix object).
            %
            
            if ~(isa(prior, 'gaussd') && prior.has_ip && prior.num == 1)
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
            
            if isfloat(isigma)
                if isscalar(isigma)
                    isigma = udmat(xd, isigma);
                else
                    if ~isequal(size(isigma), [xd xd])
                        error('gaussgm:invalidarg', ...
                            'The size of isigma is invalid.');
                    end
                    isigma = gsymat(isigma);
                end
            elseif isobject(isigma)
                if ~(isigma.d == xd && isigma.n == 1)
                    error('gaussgm:invalidarg', 'isigma is invalid.');
                end
            else
                error('gaussgm:invalidarg', 'isigma is invalid.');
            end
                
            obj.pdim = pd;
            obj.xdim = xd;
            
            obj.prior = prior;
            obj.A = A;
            obj.isigma = isigma;    
            obj.gnoise = gaussd.from_ip(0, isigma, 0, 'mp');
        end
        
        
        function LP = logpri(obj, thetas)
            % Compute the log-prior of parameters
            %
            %   LP = obj.logpri(thetas);
            %
            
            LP = obj.prior.logpdf(thetas);            
        end
            
        
        function Gs = param_models(obj, thetas)
            % Get Gaussian models with given parameters
            %
            %   Gs = obj.param_models(thetas);
            %       gets the Gaussian distribution models with given
            %       parameters
            %
            
            pd = obj.pdim;
            if ~(isfloat(thetas) && ndims(thetas) == 2 && size(thetas,1) == pd)
                error('gaussgm:loglik:invalidarg', ...
                    'params should be a pdim x m numeric matrix.');
            end
            
            A_ = obj.A;
            if isempty(A_)
                Y = thetas;
            else
                Y = A_ * thetas;
            end
                        
            Gs = gaussd.from_ip(J * Y, J);            
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
        
        
        function P = get_posterior(obj, X, w)
            % Get the posterior distribution of parameters
            %
            %   P = obj.get_posterior(X, w);
            %       
            %       computes the posterior distribution of parameters
            %       based on observed data given in X. The observed
            %       samples can be weighted by w.
            %
            %       If there are n samples, then X should be a matrix
            %       of size xdim x n. 
            %
            %       In weighted case, w should be a row vector of size
            %       1 x n.
            %
            %       In output, P is a pregaussd object.
            %
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == obj.xdim)
                error('gaussgm:invalidarg', ...
                    'X should be an d x n numeric matrix.');
            end
            
            if nargin < 3
                w = 1;
            end                

            [ha, Ja] = gcondupdate(obj.A, X, obj.isigma, w);
            P = obj.prior.inject(ha, Ja);
            
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
            %       In output, theta is a column vector of size pdim x 1,
            %
                        
            if nargin < 3
                w = 1;
            end
            theta = obj.get_posterior(X, w).get_mean();
        end                        
        
    
        function Y = pos_sample(obj, X, w, n, rstream)
            % Sampling from posterior distribution of parameters
            %
            %   Y = obj.pos_sample(X, w, n);
            %
            %       draw n samples from the posterior distribution
            %       of parameters given (weighted) observations in X.
            %       One can input w as a scalar, when all observations
            %       have the same weight.
            %
            %   Y = obj.pos_sample(X, w, n, rstream);
            %
            %       additionally specifies the random number stream
            %       from which the random numbers are generated.
            %
            
            pos = obj.get_posterior(X, w);
            
            if nargin < 5
                Y = pos.sample(n);
            else
                Y = pos.sample(n, rstream);
            end
            
        end
    end
end 

