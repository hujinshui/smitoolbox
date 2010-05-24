classdef gaussd_sampler < psampler           
    % The class for sampling from Gaussian distribution(s)
    %
    
    % Created by Dahua Lin, on Mar 23, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        mu;     % the mean vectors [d x m]
        T;      % the transformation (gtrans object)
    end
    
    methods
        
        function obj = gaussd_sampler(g, rstream)
            % constructs a sampler on a gaussian distribution with complete
            % covariance
            %
            %   obj = gauss_sampler_comp(g);
            %   obj = gauss_sampler_comp(g, rstream)
            %       constructs a sampler based on a Gaussian
            %       distribution (of class gauss_mp_comp)
            %
            %       g should be an object of class gaussd_mp.
            %
            %       rstream is the underlying random number stream,
            %       if it is omitted, the default stream will be used.
            %
            
            assert(isa(g, 'gaussd_mp'), 'gaussd_sampler:invalidarg', ...
                'g should be an object of class gaussd_sampler.');
            
            if nargin < 2
                rstream = RandStream.getDefaultStream();
            end
            
            obj.nmodels = g.nmodels;
            obj.dim = g.dim;
            obj.rstream = rstream;
            
            obj.mu = g.mu;
            obj.T = chol(g.sigma);
        end
        
        
        function X = draw(obj, i, n)
            % draw iid samples from the Gaussian distributions
            %
            %   X = obj.draw(i, n);
            %       Draws samples from selected distributions.
            %
            %       Here, i is the index array of the selected
            %       distributions, and n is the number of samples
            %       to be drawn from them.
            %       In particular, it draws n(j) samples from the
            %       distribution whose index is i(j).
            %
            %       If i is empty, then by default it is equivalent
            %       to being set to 1:n.
            %       
            %       If one would like to draw the same number of
            %       samples from each selected distribution, then
            %       n can be a scalar.
            %
            
            if isempty(i)
                i = 1 : obj.nmodels;
            else
                assert(isnumeric(i) && isvector(i), ...
                    'gaussd_sampler:draw:invalidarg', ...
                    'i should be a numeric vector.');
            end
            
            assert(isnumeric(n) && (isscalar(n) || numel(n) == numel(i)), ...
                'gaussd_sampler:draw:invalidarg', ...
                'n should be either a scalar or a numeric vector of the same size as i.');
            
            m = numel(i);
            d = obj.dim;
            
            if isscalar(n) && m > 1
                n = n * ones(1, m);
            end
            
            rs = obj.rstream;
            t = obj.T;
            mv = obj.mu;
            
            shared_cov = (obj.T.num == 1);
            
            if m == 1
                if ~shared_cov
                    X = bsxfun(@plus, t.take(i) * rs.randn(d, n), mv(:, i));
                else
                    X = bsxfun(@plus, t.take(1) * rs.randn(d, n), mv(:, i));
                end
            else
                if isscalar(n)                    
                    n = n * ones(1, m);
                end                                
                X = zeros(d, sum(n), class(mv));
                
                if ~shared_cov
                    k = 0;
                    for j = 1 : m
                        cX = bsxfun(@plus, t.take(i(j)) * rs.randn(d, n(j)), mv(:,j));
                        X(:, k+1:k+n(j)) = cX;
                        k = k + n(j);
                    end
                else
                    k = 0;
                    for j = 1 : m
                        cX = bsxfun(@plus, t.take(1) * rs.randn(d, n(j)), mv(:,j));
                        X(:, k+1:k+n(j)) = cX;
                        k = k + n(j);
                    end
                end
            end                        
        end
        
    end
end