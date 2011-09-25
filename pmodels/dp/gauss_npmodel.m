classdef gauss_npmodel < nonparam_model
    % Non-parametric Gaussian model    
    %
    
    % Created by Dahua Lin, on Sep 20, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')     
        dim;        % the space dimension
        gbase;      % the base Gaussian distribution        
        Cx;         % the covariance of conditional generation
        invCx;      % the inverse covariance
        
        gmargin;    % the marginal Gaussian distribution (w.r.t. base)
    end
    
    
    methods
                       
        function obj = gauss_npmodel(gbase, Cx)
            % Constructs a Nonparametric Gaussian model obj
            %
            %   obj = gauss_npmodel(gbase, Cx);
            %
            %       constructs a non-parametric Gaussian model object
            %       based on a Gaussian generative model
            %       
            %       - gbase:    the base Gaussian distribution
            %       - Cx:       the covariance of observation noises           
            %
            
            obj.supp_inheritance = true;
            
            if ~(isa(gbase, 'gaussd') && gbase.num == 1 && gbase.has_mp && gbase.has_ip)
                error('gauss_npmodel:invalidarg', ...
                    'gbase should be an object of gaussd with both mp & ip.');
            end
            
            if ~(is_pdmat(Cx) && Cx.n == 1)
                error('gauss_npmodel:invalidarg', ...
                    'Cx should be a pdmat struct with single matrix.');
            end
            
            if Cx.d ~= gbase.dim
                error('gauss_npmodel:invalidarg', ...
                    'The dimension of gbase and that of Cx are inconsistent.');
            end
            
            mg_cov = pdmat_plus(gbase.C, Cx);
            
            obj.dim = Cx.d;
            obj.gbase = gbase;
            obj.Cx = Cx;
            obj.invCx = pdmat_inv(Cx);
            obj.gmargin = gaussd.from_mp(gbase.mu, mg_cov, 'ip');
        end
        
        
        function n = get_num_samples(obj, X)
            % Gets the number of a given set of samples
            %
            %   n = obj.get_num_samples(X);
            %
            %       verifies the validity of X as a data set and returns
            %       the number of samples in X.
            %
            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == obj.dim)
                error('gauss_npmodel:invalidarg', ...
                    'X is not a valid sample matrix.');
            end            
            n = size(X, 2);
        end
        
        
        function a = posterior_atom(obj, X, I, pri_a)
            % Samples a new atom based on given data
            %
            %   a = obj.create_atom(X, I);
            %       samples a new atom based on the data in X, 
            %       with indices in I.
            %
            %   a = obj.create_atom(X, I, pri_a);
            %       samples a new atom based on the data selected by I,
            %       with respect to the prior given by pri_a.
            %
            %       Here, pri_a can be either a vector or a gaussd object.
            %
            
            % verify input
            
            if nargin < 4
                pri_a = obj.gbase;
            else
                d = obj.dim;
                if isnumeric(pri_a)
                    if ~(isfloat(pri_a) && isreal(pri_a) && ...
                            isequal(size(pri_a), [d 1]))
                        error('gauss_npmodel:invalidarg', 'pri_a is invalid.');
                    end
                    a = pri_a;
                    return;
                    
                elseif isa(pri_a, 'gaussd')
                    if ~(pri_a.num == 1 && pri_a.dim == d && pri_a.has_ip)
                        error('gauss_npmodel:invalidarg', 'pri_a is invalid.');
                    end
                    
                else
                    error('gauss_npmodel:invalidarg', 'pri_a is invalid.');
                end
            end
            
            % sample from posterior
            
            cX = X(:, I);
            n = size(cX, 2);
            
            Jx = obj.invCx;
            
            h = pri_a.h + pdmat_mvmul(Jx, sum(cX, 2));
            J = pdmat_plus(pri_a.J, pdmat_scale(Jx, n));
            
            C = pdmat_inv(J);
            mu = pdmat_mvmul(C, h);
            
            a = gsample(mu, C, 1);                
        end
        
        
        function L = evaluate_loglik(obj, a, X)
            % Evaluates the log-likelihood for an atom
            %
            %   L = obj.evaluate_loglik(a, X);            
            %       evaluates the log-likelihood values of all samples
            %       with respect to the input atom a.
            %
            %   L = obj.evaluate_loglik([], X);
            %   L = obj.evaluate_loglik(pri_a, X);
            %
            %       evaluates the (marginal) log-likelihood values of all 
            %       samples with respect to the base distribution or the
            %       given prior.
            %       
            
            if isempty(a)
                g = obj.gmargin;                
            elseif isa(a, 'gaussd')
                g = gaussd.from_mp(a.mu, pdmat_plus(a.C, obj.Cx), 'ip');                
            else
                g = gaussd.from_mp(a, obj.Cx, 'ip');
            end            
            L = g.logpdf(X);
        end
                
    end
    
end
