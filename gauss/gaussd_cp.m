classdef gaussd_cp < expfdistr
    % The class to represent Gaussian distribution(s) with canonical
    % parameters
    %
    % All Gaussian distribution has the following quadratic form
    %
    %   log p = <theta1, x> - 0.5 * <theta2, xx^T> - A(theta1, theta2)
    %
    %   theta1 = inv(cov) * mu
    %   theta2 = inv(cov)
    %
    %   In addition, we use
    %
    %   theta3 = mu' * inv(cov) * mu
    %
    
    % Created by Dahua Lin, on Mar 22, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
        
    properties(GetAccess='public', SetAccess='protected')       
        dim;        % the dimensionality of the underlying space
        
        theta1;     % the quantity corresponding to inv(cov) * mu
        theta2;     % the quantity corresponding to inv(cov)
        theta3;     % the quantity corresponding to mu' * inv(cov) * mu
    end
    
    methods
        function obj = gaussd_cp(t1, t2, t3)
            % constructs a Gaussian distribution object with canonical 
            % parametersa
            %
            %   obj = gaussd_cp(theta1, theta2, theta3);
            %       constructs a Gaussian distribution object with
            %       the parameters theta1, theta2, and theta3 respectively
            %
            %       If there are m distributions on a d-dimensional
            %       vector space, then 
            %       - theta1 should be a d x m matrix
            %       - theta2 should be a pdmat object with dim == d and
            %         num == m (the covariance is not shared), or 
            %         num == 1 (the covariance is shared)
            %       - theta3 should be a 1 x m row vector.
            %     
            
            if ~(isfloat(t1) && ndims(t1) == 2) 
                error('gaussd_cp:invalidarg', ...
                    'theta1 should be a numeric matrix.');
            end
            
            [d, m] = size(t1);
            
            if ~(isa(t2, 'pdmat') && t2.dim == d && (t2.num == 1 || t2.num == m))
                error('gaussd_cp:invalidarg', ...
                    'theta2 should be a pdmat object with dim == d and num being 1 or m.');
            end
            
            if ~(isfloat(t3) && isequal(size(t3), [1 m]))
                error('gaussd_cp:invalidarg', ...
                    'theta3 should be a 1 x m numeric row vector.');
            end
            
            obj.dim = d;            
            obj.nmodels = m;
            obj.logpar = 0.5 * (t3 + d * log(2*pi) - logdet(t2));
            
            obj.theta1 = t1;
            obj.theta2 = t2;
            obj.theta3 = t3;
        end
        
                
        function D = mahdist(obj, X, i)
            % Compute the Mahalanobis distances
            %
            %   D = mahdist(obj, X);
            %   D = mahdist(obj, X, i);
            %       compute the Mahalanobis distances from the samples
            %       given in X to the selected models.
            %
            %       i is a index of an array of indices of selected 
            %       models. If i is omitted, then all models are selected.
            %       
            %       If there are m selected models and n samples, then
            %       in the output, D should be a m x n matrix, with
            %       D(i, j) being the Mahalanobis distance from the 
            %       j-th sample to the i-th center w.r.t. the i-th
            %       covariance matrix.
            %
            
            if nargin < 3
                D = sqrt(max(mahdist_sq(obj, X), 0));
            else
                D = sqrt(max(mahdist_sq(obj, X, i), 0));
            end
        end
        
        function D = mahdist_sq(obj, X, i)
            % Compute the squared Mahalanobis distances
            %
            %   D = mahdist_sq(obj, X);
            %   D = mahdist_sq(obj, X, i);
            %       compute the squared Mahalanobis distances from the 
            %       samples given in X to the selected models.
            %
            %       i is a index of an array of indices of selected 
            %       models. If i is omitted, then all models are selected.
            %       
            %       If there are m selected models and n samples, then
            %       in the output, D should be a m x n matrix, with
            %       D(i, j) being the squared Mahalanobis distance from 
            %       the j-th sample to the i-th center w.r.t. the i-th
            %       covariance matrix.
            %
            
            t3 = obj.theta3;
            if nargin < 3
                Y = compute_clinterm(obj, X);
            else
                Y = compute_clinterm(obj, X, i);
                t3 = t3(i(:));
            end
            
            if isscalar(t3) || size(Y, 1) == 1
                D = (-2) * Y + t3;
            else
                D = bsxfun(@plus, (-2) * Y, t3.');
            end
        end     
        
        function LT = compute_clinterm(obj, X, i)
            % Compute the canonical linear term (of the selected models)
            %
            %   LT = compute_clinterm(obj, X);
            %   LT = compute_clinterm(obj, X, i);
            %       compute the canonical linear term of selectd models.
            %
            %       i is a index of an array of indices of selected 
            %       models. If i is omitted, then all models are selected.
            %
            %       If there are m models selected and n samples in X,
            %       then LT is a m x n matrix, with LT(i, j) corresponding 
            %       to i-th model and j-th sample.
            %
            
            % verify input arguments 
            d = obj.dim;
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == d)
                error('gauss_cp:compute_clinterm:invalidarg', ...
                    'X should be a numeric matrix of size d x n.');
            end
            
            % compute
            if nargin < 3
                i = [];
            end
            
            % linear part
            
            if isempty(i)
                LT1 = obj.theta1' * X;
            else
                LT1 = obj.theta1(:,i)' * X;
            end
            
            % quadratic part
            
            if isempty(i) || obj.theta2.num == 1
                LT2 = obj.theta2.quadterm(X);
            else
                LT2 = obj.theta2.take(i).quadterm(X);
            end
                
            if size(LT2, 1) == 1 && size(LT1, 1) > 1
                LT = bsxfun(@minus, LT1, 0.5 * LT2);
            else
                LT = LT1 - 0.5 * LT2;
            end                                    
        end
        
                
        function L = compute_logbase(obj, X) %#ok<INUSD,MANU>
            % Evaluate the canonical linear term on given samples
            L = 0;
        end        
        
        
        function R = to_mp(obj)
            % Converts the object to mean-parameterized object
            %
            %   R = to_mp(obj);
            %
            
            sigma = inv(obj.theta2);
            if sigma.num == 1
                mu = sigma * obj.theta1; %#ok<MINV>
            else
                mu = cmult(sigma, obj.theta1);
            end
            
            R = gaussd_mp(mu, sigma);
        end
        
        
        function R = conj_update(obj, U1, U2)
            % Performs conjugate update of the model
            %
            %   R = conj_update(obj, U1, U2);
            %       updates the model in obj by adding U1 to theta1
            %       and adding U2 to theta2.
            %
            %       If obj.num == 1, then U1 and U2 can represent
            %       any number of models, say m. Then R is also
            %       a gaussd_cp object, which contains m models.
            %
            %       If obj.num == m > 1, then U1 and U2 should 
            %       represent m models. The update are done 
            %       correspondingly.
            %
            
            % verify input
            
            d = obj.dim;
            m0 = obj.num;
            
            if ~(isfloat(U1) && isreal(U1) && ndims(U1) == 2 && size(U1, 1) == d)
                error('gaussd_cp:conj_update:invalidarg', ...
                    'U1 should be a real matrix with size(U1,1) == obj.dim.');
            end
            
            if ~(isa(U2, 'pdmat') && is_subtype(obj.theta2, U2) && U2.dim == d)
                error('gaussd_cp:conj_update:invalidarg', ...
                    'U2 should be a pdmat object compatible with theta2 with dim d.');
            end
            
            if size(U1, 2) ~= U2.num
                error('gaussd_cp:conj_update:invalidarg', ...
                    'U1 and U2 represent different numbers of models.');
            end
            m1 = size(U1, 2);
            
            if ~(m0 == 1 || m0 == m1)
                error('gaussd_cp:conj_update:invalidarg', ...
                    'The number of models in object does not match that for U1 and U2.');
            end
            
            % compute
            
            if m0 == m1
                T1 = obj.theta1 + U1;
            else
                T1 = bsxfun(@plus, obj.theta1, U1);
            end
            
            T2 = obj.theta2 + U2;
                    
            mu = cldiv(T2, T1);
            R = gaussd_cp(T1, T2, sum(mu .* T1, 1));
        end
    end
      
    
    methods(Static)
        function obj = from_mean_and_cov(mu, sigma)
            % Constructs a canonically parameterized Gaussian object
            % from the given mean and covariance
            %
            %   obj = gaussd_cp.from_mean_and_cov(mu, sigma);
            %       constructs a Gaussian object with canonical 
            %       parameterization using the mean given in mu
            %       and the covariance given by sigma.
            %   
            %       If there are m distributions in a d-dimensional
            %       space, then mu should be a d x m matrix, and
            %       sigma should be a pdmat object with dim == d
            %       and num == m (non-shared covariance) or 
            %       num == 1 (shared covariance).
            %
            
            if ~(isfloat(mu) && ndims(mu) == 2)
                error('gaussd_cp:from_mean_and_cov:invalidarg', ...
                    'mu should be a numeric matrix.');
            end
            
            [d, m] = size(mu);
            
            if ~(isa(sigma, 'pdmat') && sigma.dim == d && ...
                (sigma.num == 1 || sigma.num == m))
                error('gaussd_cp:from_mean_and_cov:invalidarg', ...
                    'sigma should be a pdmat object with dim == d and num being 1 or m.');
            end
            
            t2 = inv(sigma);
            if t2.num == 1
                t1 = t2 * mu; %#ok<MINV>
            else
                t1 = cmult(t2, mu);                            
            end
            t3 = sum(t1 .* mu, 1);
            
            obj = gaussd_cp(t1, t2, t3);
        end
        
        function obj = from_mean_and_icov(mu, isigma)
            % Constructs a canonically parameterized Gaussian object
            % from the given mean and covariance
            %
            %   obj = gaussd_cp.from_mean_and_icov(mu, isigma);
            %       constructs a Gaussian object with canonical 
            %       parameterization using the mean given in mu
            %       and the inverse covariance given by isigma.
            %
            
            if ~(isfloat(mu) && ndims(mu) == 2)
                error('gaussd_cp:from_mean_and_cov:invalidarg', ...
                    'mu should be a numeric matrix.');
            end
            
            [d, m] = size(mu);
            
            if ~(isa(isigma, 'pdmat') && isigma.dim == d && ...
                (isigma.num == 1 || isigma.num == m))
                error('gaussd_cp:from_mean_and_cov:invalidarg', ...
                    'sigma should be a pdmat object with dim == d and num being 1 or m.');
            end
            
            t2 = isigma;
            if t2.num == 1
                t1 = t2 * mu; 
            else
                t1 = cmult(t2, mu);                            
            end
            t3 = sum(t1 .* mu, 1);
            
            obj = gaussd_cp(t1, t2, t3);
        end
    end
end

