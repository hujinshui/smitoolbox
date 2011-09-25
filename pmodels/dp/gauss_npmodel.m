classdef gauss_npmodel < nonparam_model
    % Non-parametric Gaussian model
    %
    % This is a wrapper of Gaussian generative model to nonparam_model
    % interfaces.
    %
    
    properties(GetAccess='public', SetAccess='private')
        gm;         % the Gaussian generative model
        pri_mu;     % the prior mean vector
        
        gmargin;    % the marginal Gaussian model
        
        X;          % the observed data
    end
    
    
    methods
       
        function obj = gauss_npmodel(gm, pri_mu, X)
            % Constructs a Nonparametric Gaussian model obj
            %
            %   obj = gauss_npmodel(gm, pri_mu, X);
            %
            %       constructs a non-parametric Gaussian model object
            %       based on a Gaussian generative model
            %       
            %       - gm:   the underlying Gaussian generative model
            %       - X:    the observed data, size d x n
            %
            
            if ~(isa(gm, 'gaussgm') && ~isempty(gm.Cx) && ~isempty(gm.Cu))
                error('gauss_npmodel:invalidarg', ...
                    'gm should be a gaussgm object with non-empty Cx and Cu.');
            end
            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
                error('gauss_npmodel:invalidarg', ...
                    'X should be a real matrix of size d x n.');
            end
            
            if ~(isfloat(pri_mu) && isreal(pri_mu) && ...
                    isequal(size(pri_mu), [gm.udim, 1]))
                error('gauss_npmodel:invalidarg', ...
                    'pri_mu should be a real vector of size udim x 1.');
            end
            
            if size(X, 1) ~= gm.xdim
                error('gauss_npmodel:invalidarg', ...
                    'X should have size(X,1) == gm.xdim.');
            end
            
            obj = obj@nonparam_model(size(X, 2));
            obj.supp_inheritance = false;
            obj.gm = gm;
            obj.pri_mu = pri_mu;
            
            mg_cov = pdmat_plus(gm.Cu, gm.Cx);
            obj.gmargin = gaussd.from_mp(pri_mu, mg_cov, 'ip');
            
            obj.X = X;            
        end
        
        
        function a = posterior_atom(obj, I)
            % Samples a new atom based on given data
            %
            %   a = obj.create_atom(I);
            %       samples a new atom based on the data with indices 
            %       in I.
            %
            
            a = obj.gm.sample_u(obj.X(:, I), [], obj.pri_mu);
        end
        
        
        function L = evaluate_loglik(obj, a)
            % Evaluates the log-likelihood for an atom
            %
            %   L = obj.evaluate_loglik(a);
            %       evaluates the log-likelihood values of all samples
            %       with respect to the input atom a.
            %
            %   L = obj.evaluate_loglik();
            %       evaluates the (marginal) log-likelihood values of all 
            %       samples with respect to the base distribution
            %       
            
            if nargin < 2 || isempty(a)
                L = obj.gmargin.logpdf(obj.X);
            else
                L = obj.gm.loglik(a, obj.X);
            end                        
        end
                
    end
    
end
