classdef wishartd
    % The class to represent Wishart distribution
    %
    %   Wishart distribution is a multi-dimensional generalization of
    %   gamma distribution. 
    %
    %   A Wishart distribution is characterized by two parameters:
    %   - deg:      degree of freedom
    %   - S:        scale matrix
    %
    %   Currently, the class only supports single-distribution objects,
    %   which can be one or multiple dimensional.
    %
    
    % Created by Dahua Lin, on Sep 1, 2011
    %
    
    
    %% Properties
    
    properties
        num = 1;    % the number of models contained in the object
        dim;        % dimension (sample matrix size: dim x dim)
        
        deg;        % degree of freedom ( > p-1 )
        S;          % the scale matrix  (pdmat struct)
        
        lpconst;    % pre-computed const term for logpdf evaluation
                    % 
                    % = q*d*log(2) + q * log|S| + log(Gamma_d(q))
                    %   
                    %   with q = deg / 2
                    
        invS;       % the inverse of S.
    end
    
    
    %% Construction
    
    methods
        
        function obj = wishartd(S, deg, op)
            % Constructs a Wishart distribution object
            %
            %   obj = wishartd(S, deg);
            %       constructs a Wishart distribution of given degree
            %       of freedom (deg) and scale matrix (S).
            %
            %   obj = wishartd(S, deg, 'pre');
            %
            %       performs the construction and pre-computes lpconst
            %       and invS to accelerate evaluation of logpdf or pdf.
            %
              
            % verify inputs
            
            if ~(is_pdmat(S) && S.n == 1)
                error('wishartd:invalidarg', ...
                    'S should be a pdmat struct with S.n == 1.');
            end
            
            if ~(isnumeric(deg) && isscalar(deg) && isreal(deg) && deg >= 0)
                error('wishartd:invalidarg', ...
                    'deg should be a non-negative real value.');
            end
            deg = double(deg);
                                    
            pre_comp = 0;
            if nargin >= 3
                if ~(ischar(op) && strcmpi(op, 'pre'))
                    error('wishartd:invalidarg', ...
                        'The 3rd arg can only be ''pre''.');
                end
                pre_comp = 1;
            end
            
            % create the object
            
            obj.dim = S.d;
            obj.deg = deg;
            obj.S = S;
                        
            if pre_comp
                obj.lpconst = wishartd.calc_lpconst(deg, S);
                obj.invS = pdmat_inv(S);
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
            
            m = obj.deg;
            if m == 1
                M = obj.S;
            else
                M = pdmat_scale(obj.S, m);
            end
        end
        
    
        function M = mode(obj)
            % Gets the mode of the distribution
            %
            %   M = mode(obj);
            %       returns the mode matrix in form of a pdmat struct
            %
            
            m = obj.deg;
            d = obj.dim;
            
            M = pdmat_scale(obj.S, max(m-d-1, 0));
        end
        
    end
    
    
    %% Evaluation
    
    methods
        
        function L = logpdf(obj, W)
            % Evaluates the log-pdf 
            %
            %   L = obj.logpdf(W);
            %       evaluates the logarithm of pdf of sample matrix(ices)
            %       given in form of a pdmat struct.
            %   
            
            d = obj.dim;
            if ~(is_pdmat(W) && W.d == d)
                error('wshartd:invalidarg', ...
                    'W should be a pdmat struct with W.d == d.');
            end
            
            m = obj.deg;
            
            lpc = obj.lpconst;
            if isempty(lpc)
                lpc = wishartd.calc_lpconst(m, obj.S);
            end
            
            J = obj.invS;
            if isempty(J)
                J = pdmat_inv(obj.S);
            end
            
            t1 = (0.5 * (m - d - 1)) * pdmat_lndet(W);
            t2 = 0.5 * pdmat_dot(J, W);
            
            L = t1 - t2 - lpc;
        end
        
        
        function P = pdf(obj, W)
            % Evaluates the pdf values
            %
            %   P = obj.pdf(W);
            %       evaluates the pdf values of the matrices in W,
            %       which is a pdmat struct.
            %
            
            P = exp(logpdf(obj, W));            
        end
    
    end
    
    
    methods(Static, Access='private')
        
        function lpc = calc_lpconst(deg, S)
            % compute the logpdf const term
            
            q = deg * 0.5;
            ln2 = 0.69314718055994530941723212;
            d = S.d;
            
            lpc = q * (d * ln2 + pdmat_lndet(S)) + mvgammaln(d, q);
        end                
    end    
    
    
    %% Sampling
    
    methods
        
        function W = sample(obj, n)
            % Samples from the distribution
            %
            %   W = obj.sample();
            %   W = obj.sample(n);
            %
            %       Draws n (if omitted, n = 1) samples from the 
            %       distribution, and packs them in a pdmat struct W.
            %

            if nargin < 2
                n = 1;
            end
            
            m = obj.deg;
            S_ = obj.S;
            
            if isequal(S_.v, 1)
                W = wishart_sample(S_.d, m, n);
            else
                W = wishart_sample(S_, m, n);
            end
            
            W = pdmat('f', S_.d, W);
        end
    
    end
    
    
end



