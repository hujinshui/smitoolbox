classdef scale_invchi2d
    % The class to represent Scaled-Inverse Chi-Square distribution
    %
    %   A scaled inverse chi-square distribution is essentially a
    %   re-parametrization of inverse gamma distribution, that is
    %   more convenient in serving as the conjugate prior of normal
    %   variance.
    %
    %   The distribution is characterized by two parameters:
    %
    %   - nu:       degree of freedom
    %   - sigma2:   inverse scale
    %
    %   The relation between this and inverse gamma distribution is:
    %
    %     Scaled-Inv-Chi-Square(nu, sigma2) ~ Inv-Gamma(alpha, beta).
    %
    %   with 
    %
    %     alpha = nu / 2;
    %     beta  = nu * sigma2 / 2;
    %
    %   Currently, this class only supports single-distribution objects,
    %   which, however, can be one or multi-dimensional.
    %
    
    % Created by Dahua Lin, on Sep 1, 2011
    %
    
    %% Properties
    
    properties
        dim;
        num = 1;
        
        nu;         % the degree of freedom
        sigma2;     % base variance (inverse scale)
        
        alpha;      % inverse-gamma shape = nu / 2
        beta;       % inverse-gamma scale = nu * sigma2 / 2;
        lpconst;    % the constant in logpdf evaluation
    end
    
    %% Construction
    
    methods
        
        function obj = scale_invchi2d(nu, sigma2, d)
            % Constructs a scaled inverse chi-square distribution object
            %
            %   obj = scale_invchi2d(nu, sigma2);
            %   obj = scale_invchi2d(nu, sigma2, d);
            %
            %       constructs a scaled inverse chi-square distribution
            %       with nu deg. of freedom and base variance sigma2.
            %
            %       For a d-dimensional distribution,
            %       - nu:       a scalar
            %       - sigma2:   a scalar or a d x 1 column vector
            %
            %       Note when sigma2 is a scalar, and d > 1, then
            %       the 3rd argment must be explicitly given,
            %       otherwise, it will consider d == 1.
            %
            
            if ~(isfloat(nu) && isreal(nu) && isscalar(nu) && nu > 0)
                error('scale_invchi2d:invalidarg', ...
                    'nu should be a numeric scalar.');
            end
            
            if ~(isfloat(sigma2) && isreal(sigma2) && ...
                    ndims(sigma2) == 2 && size(sigma2,2) ==1)
                error('scale_invchi2d:invalidarg', ...
                    'sigma2 should be either a scalar or a column vector.');
            end
            d_ = size(sigma2, 1);
            
            if nargin < 3
                d = d_;
            else
                if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1)
                    error('scale_invchi2d:invalidarg', ...
                        'd should be a positive integer scalar.');
                end
                if ~(d_ == 1 || d_ == d)
                    error('scale_invchi2d:invalidarg', ...
                        'Inconsistency between sigma2 and d.');
                end
            end
            
            % create object
            
            a = 0.5 * nu;
            b = a * sigma2;
            
            obj.dim = d;
            obj.nu = nu;
            obj.sigma2 = sigma2;
            
            obj.alpha = a;
            obj.beta = b;
            
            if d == 1
                lpc = a * log(b) - gammaln(a);
            else
                if d_ == 1
                    lpc = (a * log(b) - gammaln(a)) * d;
                else
                    lpc = a * sum(log(b), 1) - gammaln(a) * d;
                end
            end            
            obj.lpconst = lpc;
        end
        
    end
    
    %% Model statistics
    
    methods
        
        function v = mean(obj)
            % Evaluates the mean(s) of the distribution
            %
            %   v = mean(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding
            %       to the i-th distribution in obj.
            %
            %   Note this function applies only to the case with
            %   alpha > 1.
            %
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            v = b ./ (a - 1);
            if d > size(b, 1)
                v = v(ones(d, 1), :);
            end
        end
        
        function v = var(obj)
            % Evaluates the variance(s) of the distribution
            %
            %   v = var(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding
            %       to the i-th distribution in obj.
            %
            %   Note this function applies only to the case with
            %   alpha > 2.
            %
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            v = (b.^2) ./ ((a - 1).^2 .* (a - 2));
            if d > size(b, 1)
                v = v(ones(d, 1), :);
            end
        end
        
        function v = mode(obj)
            % Gets the mode of the distribution
            %
            %   v = mode(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding
            %       to the i-th distribution in obj.
            %
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            v = b ./ (a + 1);
            if d > size(b, 1)
                v = v(ones(d, 1), :);
            end
        end
        
        
        function v = entropy(obj)
            % Evaluates the entropy value(s) of the distribution
            %
            %   v = entropy(obj);
            %
            %       v will be a 1 x n matrix, with v(:,i) corresponding
            %       to the i-th distribution in obj.
            %
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            v = (a + gammaln(a) - (1 + a) .* psi(a)) + log(b);
            if d > 1
                if size(b, 1) == 1
                    v = d * v;
                else
                    v = sum(v, 1);
                end
            end
 
        end
        
    end
    
    
    %% Evaluation
    
    methods
    
        function L = logpdf(obj, X)
            % Evaluates the logarithm of pdf 
            %
            %   L = obj.logpdf(X);
            %
            %       evaluates the logarithm of pdf at the samples given
            %       as columns of X.
            %
            
            d = obj.dim;
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == d)
                error('scale_invchi2d:invalidarg', ...
                    'X should be a numeric matrix with size(X,1) == dim.');
            end
            
            a = obj.alpha;
            b = obj.beta;
            lpc = obj.lpconst;
            
            if d == 1
                L = -( (a + 1) * log(X) + b ./ X) + lpc;
            else
                if size(b, 1) == 1
                    L = -((a + 1) * sum(log(X), 1) + b * sum(1./X, 1)) + lpc;
                else
                    L = -((a + 1) * sum(log(X), 1) + b' * (1 ./ X)) + lpc;                    
                end
            end
        end
        
        
        function P = pdf(obj, X)
            % Evaluates the pdf values
            %
            %   P = obj.pdf(X);
            %
            %       evaluates the pdf values at columns of X.
            %
            
            P = exp(logpdf(obj, X));            
        end
    
    end
    
    
    %% Sampling
    
    methods
        
        function X = sample(obj, n)
            % Samples from the distribution
            %
            %   X = obj.sample(n);
            %
            %       Draws n samples from the scaled inverse chi-square
            %       distribution.
            %
            
            X = gamma_sample(obj.alpha, 1./obj.beta, n, obj.dim);    
            X = 1 ./ X;
        end
    
    end
    
end
