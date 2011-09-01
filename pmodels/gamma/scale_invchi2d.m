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
            obj.lpconst = a * log(b) - gammaln(a);            
        end
    
    end
    
    
end
