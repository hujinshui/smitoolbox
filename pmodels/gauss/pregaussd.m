classdef pregaussd
    % The class to capture 1st and 2nd order information params of Gaussian
    %
    % This is a light-weight class basically for efficient computation
    % related to posterior Gaussian.
    %
    
    % Created by Dahua Lin, on Sep 16, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        dim;    % the space dimension
        num;    % the number of models contained in the object
        h;      % the potential vector
        J;      % the information matrix
    end
    
    methods
        
        function G = pregaussd(h, J)
            % Construct a pregaussd object
            %
            %   G = pregaussd(h, J);
            %
            
            if ~(isfloat(h) && ndims(h) == 2)
                error('pregaussd:invalidarg', 'h should be a numeric matrix.');
            end
            
            [d, n] = size(h);
            
            if ~(isobject(J) && J.d == d && (J.n == n || J.n == 1))
                error('pregaussd:invalidarg', 'J is invalid.');
            end
            
            G.dim = d;
            G.num = n;
            G.h = h;
            G.J = J;
        end
        
        
        function G1 = to_gaussd(G, ump)
            % Make a gaussd object from this
            %
            %   G1 = G.to_gaussd();
            %   G1 = G.to_gaussd('mp');
            %
            
            if nargin < 2
                G1 = gaussd.from_ip(G.h, G.J);
            else
                G1 = gaussd.from_ip(G.h, G.J, [], ump);
            end                
        end
        
        
        function v = get_mean(G)
            % compute the mean vector(s) of the Gaussian distribution
            %
            %   v = G.get_mean();
            %   
            
            v = G.J \ G.h;
        end
        
        
        function C = get_cov(G)
            % compute the covariance object of the Gaussian distribution
            %
            %   C = G.get_cov();
            %
            
            C = inv(G.J);
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
                                    
            h_ = G.h;
            J_ = G.J;
                                    
            if nargin < 3
                X = randn(G.dim, n);
            else
                X = rstream.randn(G.dim, n);
            end
            
            C = inv(J_);
            X = C.choltrans(X);
            
            if ~isequal(h_, 0)
                mu = C * h_; %#ok<MINV>
                X = bsxfun(@plus, X, mu);
            end
        end
            
    end
end

