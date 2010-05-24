classdef gcondmap
    % The class for MAP estimation of Gaussian conditional models
    %
    
    % Created by Dahua Lin, on Apr 13, 2010
    %
    
    properties(Constant)
        est_type = 'map';
    end
    
    properties(GetAccess='public', SetAccess='private')
        dim;            % dimension of vector space of the variables being conditioned
        gpri;           % Gaussian prior 
        gcond;          % Gaussian conditional model
    end
    
    methods
        
        function obj = gcondmap(gpri, gcond)
            % Constructs a Gaussian MAP estimator
            %
            %   obj = gaussd_map(gpri, gcond);
            %       constructs an MAP estimator for Gaussian conditional
            %       model. The prior model and conditional likelihood
            %       model are respectively given by gpri and gcond.
            %
            
            assert((isa(gpri, 'gaussd_cp') || isa(gpri, 'gaussd_mp')) && ...
                gpri.nmodels == 1, ...
                'gcondmap:invalidarg', ...
                'gpri should be a Gaussian model object with single model.');
            
            if isa(gpri, 'gaussd_mp')
                gpri = to_cp(gpri);
            end
            
            assert(isa(gcond, 'gcondmodel') && gcond.cdim == gpri.dim, ...
                'gcondmap:invalidarg', ...
                'gcond should be an object of class gcondmodel with cdim == d.');            
            
            obj.dim = gpri.dim;
            obj.gpri = gpri;
            obj.gcond = gcond;                    
        end
        
        function Z = estimate(obj, X, W)
            % Estimate the Gaussian mean in MAP manner
            %
            %   Z = obj.estimate(X);
            %   Z = obj.estimate(X, W);
            %       performs maximum a posteriori estimation of hidden 
            %       variables with fixed covariance.
            %
                                    
            g0 = obj.gpri;
            gc = obj.gcond;
            
            [U1, U2] = gc.conj_update(X, W);
                        
            if size(U1,2) == 1
                T1 = g0.theta1 + U1;
            else
                T1 = bsxfun(@plus, g0.theta1, U1);
            end
            T2 = g0.theta2 + U2;
            
            Z = cldiv(T2, T1);                                 
        end
        
        
        function lpri = eval_logpri(obj, Z)
            % Evaluate log-prior of the input models
            %
            %   lpri = eval_logpri(obj, Z);
            %       Evaluates the logarithm of prior pdf value of the
            %       hidden variables.
            %
            %       In the input, Z is a matrix of size cdim x m, and then
            %       lpri is of size 1 x m.
            %
            
            lpri = logprob(obj.gpri, Z);         
        end
        
        
        function L = eval_loglik(obj, Z, X)
            % Evaluate the log likelihood of samples for given models
            %
            %   L = eval_loglik(obj, Z, X);
            %       evaluates the log likelihood of the samples in X
            %       with respect to the given hidden variables Z.
            %
            %       Suppose the size of Z is cdim x m, and then L is of 
            %       size m x n.
            %
            
            L = eval_loglik(obj.gcond, Z, X);           
        end
        
    end
    
end