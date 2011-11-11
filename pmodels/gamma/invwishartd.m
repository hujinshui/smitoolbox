classdef invwishartd
    % The class to represent inverse Wishart distribution
    %
    %   Inverse Wishart distribution is a multi-dimensional generalization 
    %   of inverse gamma distribution. 
    %
    %   An inverse Wishart distribution is characterized by two parameters:
    %   - deg:      degree of freedom
    %   - S:        inverse scale matrix
    %
    %   Currently, the class only supports single-distribution objects,
    %   which can be one or multiple dimensional.
    %
    %   The relation between Inverse Wishart and Wishart is
    %
    %   C ~ InvWishart(Phi, m) <--> inv(C) ~ Wishart(inv(Phi), m)
    %
    %   Inverse Wishart distribution is useful in Bayesian modeling as
    %   a conjugate prior of Gaussian covariance matrix.
    %
    
    % Created by Dahua Lin, on Sep 2, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')
        
        num = 1;    % the number of distributions in object
        dim;        % the dimension (matrix size dim x dim)
        
        deg;        % the degree of freedom
        Phi;        % the inverse scale matrix
        
        lpconst;    % the constant for computing loglik
                    %
                    % = q*d*log(2) - q * log|Phi| + log(Gamma_d(q))
                    % 
                    % with q = deg / 2
                    %
                    
        invPhi;     % inverse of Phi, which is useful for sampling        
    end
    
    
    %% Construction
    
    methods
        
        function obj = invwishartd(Phi, deg, op)
            % Construct an inverse Wishart distribution object
            %
            %   obj = invwishartd(Phi, deg);
            %
            %       constructs an inverse Wishart distribution with
            %       inverse scale matrix Phi and the specified degree
            %       of freedom.
            %
            %       Phi should be input as a pdmat struct.
            %
            %   obj = invwishartd(Phi, deg, 'pre');
            %
            %       performs pre-computation in construction, which 
            %       will accelerate the evaluation of logpdf and sampling.
            %       
            
            % verify inputs
            
            if ~(is_pdmat(Phi) && Phi.n == 1)
                error('invwishartd:invalidarg', ...
                    'S should be a pdmat struct with S.n == 1.');
            end
            
            if ~(isnumeric(deg) && isscalar(deg) && isreal(deg) && deg >= 0)
                error('invwishartd:invalidarg', ...
                    'deg should be a non-negative real value.');
            end
            deg = double(deg);
                                    
            pre_comp = 0;
            if nargin >= 3
                if ~(ischar(op) && strcmpi(op, 'pre'))
                    error('invwishartd:invalidarg', ...
                        'The 3rd arg can only be ''pre''.');
                end
                pre_comp = 1;
            end
            
            % create the object
            
            obj.dim = Phi.d;
            obj.deg = deg;
            obj.Phi = Phi;
            
            if pre_comp
                obj.lpconst = invwishartd.calc_lpconst(deg, Phi);
                obj.invPhi = pdmat_inv(Phi);
            end            
        end                        
    
    end
    
    
    %% Statistics
    
    methods
    
        function M = mean(obj)
            % Gets the mean matrix of the distribution
            %
            %   M = mean(obj);
            %       returns the mean matrix in form of a pdmat struct.
            %
            
            s = obj.deg - obj.dim - 1;
            if s == 1
                M = obj.Phi;
            else
                M = pdmat_scale(obj.Phi, 1/s);
            end
        end
        
    
        function M = mode(obj)
            % Gets the mode of the distribution
            %
            %   M = mode(obj);
            %       returns the mode matrix in form of a pdmat struct
            %
            
            s = obj.deg + obj.dim + 1;
            M = pdmat_scale(obj.Phi, 1/s);
        end
        
    end
    
    
    %% Evaluation
    
    methods
        
        function L = logpdf(obj, C)
            % Evaluates the log-pdf 
            %
            %   L = obj.logpdf(C);
            %       evaluates the logarithm of pdf of sample matrix(ices)
            %       given in form of a pdmat struct.
            % 
            
            d = obj.dim;
            m = obj.deg;
            
            if ~(is_pdmat(C) && C.d == d)
                error('invwshartd:invalidarg', ...
                    'C should be a pdmat struct with W.d == d.');
            end
                                    
            lpc = obj.lpconst;
            if isempty(lpc)
                lpc = invwishartd.calc_lpconst(m, obj.Phi);
            end            
            
            t1 = (m+d+1) * pdmat_lndet(C);
            t2 = pdmat_dot(pdmat_inv(C), obj.Phi);
            
            L = -0.5 * (t1 + t2) - lpc;                        
        end
        
        function P = pdf(obj, C)
            % Evaluates the pdf values
            %
            %   P = obj.pdf(C);
            %       evaluates the pdf values of the matrices in C,
            %       which is a pdmat struct.
            %
            
            P = exp(logpdf(obj, C));             
        end
    end
    
    
    methods(Static, Access='private')
        
        function lpc = calc_lpconst(deg, Phi)
            % compute the logpdf const term
            
            q = deg * 0.5;
            ln2 = 0.69314718055994530941723212;
            d = Phi.d;
            
            lpc = q * (d * ln2 - pdmat_lndet(Phi)) + mvgammaln(d, q);
        end                
    end    
    
    
    %% Sampling
    
    methods
        
        function C = sample(obj, n)
            % Samples from the distribution
            %
            %   C = obj.sample();
            %   C = obj.sample(n);
            %
            %       Draws n (if omitted, n = 1) samples from the 
            %       distribution, and packs them in a pdmat struct C.
            %

            if nargin < 2
                n = 1;
            end
            
            d = obj.dim;
            m = obj.deg;
            Phi_ = obj.Phi;
            invPhi_ = obj.invPhi;
            
            % sample from Wishart
            
            if isequal(Phi_.v, 1)
                W = wishart_sample(d, m, n);
            else
                if isempty(invPhi_)
                    invPhi_ = pdmat_inv(Phi_);
                end
                W = wishart_sample(invPhi_, m, n);
            end
            
            % inverse the samples
            
            if d == 1
                C = 1 ./ W;
            elseif d == 2
                C = inv2x2(W);
            else
                if n == 1
                    C = inv(W);
                else
                    C = zeros(d, d, n, class(W));
                    for i = 1 : n
                        C(:,:,i) = inv(W(:,:,i));
                    end
                end
            end
            
            % create pdmat struct
                
            C = pdmat('f', d, C);
        end
    
    end
    
    
end


