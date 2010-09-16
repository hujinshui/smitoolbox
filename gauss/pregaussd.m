classdef pregaussd
    % The class to capture 1st and 2nd order canonical params of Gaussian
    %
    % This is a light-weight class basically for efficient computation
    % related to posterior Gaussian.
    %
    
    % Created by Dahua Lin, on Sep 16, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        dim;    % the space dimension
        num;    % the number of models contained in the object
        coef1;  % the coefficient vector of the first order term
        coef2;  % the coefficient matrix of the second order term
    end
    
    methods
        
        function G = pregaussd(c1, c2)
            % Construct a pregaussd object
            %
            %   G = pregaussd(c1, c2);
            %
            
            if ~(isfloat(c1) && ndims(c1) == 2)
                error('pregaussd:invalidarg', 'c1 should be a numeric matrix.');
            end
            
            [d, n] = size(c1);
            
            if ~(isobject(c2) && c2.d == d && (c2.n == n || c2.n == 1))
                error('pregaussd:invalidarg', 'c2 is invalid.');
            end
            
            G.dim = d;
            G.num = n;
            G.coef1 = c1;
            G.coef2 = c2;
        end
        
        
        function G1 = to_gaussd(G, ump)
            % Make a gaussd object from this
            %
            %   G1 = G.to_gaussd();
            %   G1 = G.to_gaussd('mp');
            %
            
            if nargin < 2
                G1 = gaussd.from_cp(G.coef1, G.coef2);
            else
                G1 = gaussd.from_cp(G.coef1, G.coef2, [], ump);
            end                
        end
        
        
        function v = get_mean(G)
            % compute the mean vector(s) of the Gaussian distribution
            %
            %   v = G.get_mean();
            %   
            
            v = G.coef2 \ G.coef1;
        end
        
        
        function C = get_cov(G)
            % compute the covariance object of the Gaussian distribution
            %
            %   C = G.get_cov();
            %
            
            C = inv(G.coef2);
        end
        
        
        function X = sample(G, n, rstream)
            % Samples from the Gaussian distribution
            %
            %   X = G.sample(n);
            %   X = G.sample(n, rstream);
            %
            %   Note that this function only works when num == 1
            %
            
            if G.num ~= 1
                error('pregaussd:invalidarg', 'G should contain only one distribution.');
            end
                                    
            c1 = G.coef1;
            c2 = G.coef2;
                                    
            if nargin < 3
                X = randn(G.dim, n);
            else
                X = rstream.randn(G.dim, n);
            end
            
            C = inv(c2);
            X = C.choltrans(X);
            
            if ~isequal(c1, 0)
                mu = C * c1; %#ok<MINV>
                X = bsxfun(@plus, X, mu);
            end
        end
            
    end
end

