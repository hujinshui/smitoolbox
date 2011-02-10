classdef sgmrf
    % The class to represent a Gaussian MRF with each node being a scalar
    %    
    %       
    
    % Created by Dahua Lin, on Nov 2, 2010
    % Modified by Dahua Lin, on Feb 10, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')
        nv;     % the number of nodes
        ne;     % the number of edges
        
        Jdv;    % the diagonal entries of J [nv x 1]
        Jgr;    % the information graph (capturing off-diagonal entries)
        
        es;     % the source vertices [ne x 1]
        et;     % the target vertices [ne x 1]
        ew;     % the edge weights [ne x 1]
    end
        
    methods
        
        function obj = sgmrf(v1, v2)
            % Construct a Gaussian MRF
            %
            %   obj = sgmrf(J);
            %       constructs a Gaussian MRF from the information matrix.
            %       
            %   obj = sgmrf(jdv, jg);
            %       constructs a Gaussian MRF from the diagonal values
            %       in jdv, and an off-diagonal graph jg.
            %   
            
            if nargin == 1
                J = v1;                
                if ~(isfloat(J) && ndims(J) == 2 && size(J,1) == size(J,2))
                    error('sgmrf:invalidarg', 'J should be a symmetric matrix.');
                end
                
                jdv = full(diag(J));
                [s, t, w] = find(J);
                i = find(s < t);
                jgr = gr_adjlist.from_edges('u', size(J,1), s(i), t(i), w(i));
                
            else
                jdv = v1;
                jgr = v2;
                
                if ~(isa(jgr, 'gr_adjlist') && jgr.dtype == 'u')
                    error('sgmrf:invalidarg', 'jgr should be an undirected graph.');
                end
                
                if ~(isfloat(jdv) && isequal(size(jdv), [jgr.nv, 1]))
                    error('sgmrf:invalidarg', 'jdv should be a vector of size nv x 1.');
                end
            end
            
            obj.nv = jgr.nv;
            obj.ne = jgr.ne;
            
            obj.Jdv = jdv;
            obj.Jgr = jgr;    
            
            m = jgr.ne;
            obj.es = double(jgr.es(1:m) + 1);
            obj.et = double(jgr.et(1:m) + 1);
            obj.ew = jgr.ew(1:m);
        end
        
        
        function J = infomat(obj)
            % Get the information matrix J
            %
            %   J = infomat(obj);
            %
            % Note the returned matrix J is an nv x nv sparse matrix.
            %
                                    
            n = obj.nv;
            s = obj.es;
            t = obj.et;
            w = obj.ew;            
            
            i = [(1:n)'; s; t];
            j = [(1:n)'; t; s];
            v = [obj.Jdv; w; w];
            if ~isa(v, 'double')
                v = double(v);
            end
            
            J = sparse(i, j, v, n, n);
        end
        
                             
        function v = qenergy(obj, s2, sst, mu)
            % Compute quadratic energy using covariance values
            %
            %   v = qenergy(obj, s2, sst);
            %
            %       This syntax returns v = (1/2) * tr(J * C).
            %
            %       Here, in the input arguments, s2 is the marginal
            %       variance vector and sst is the covariance values
            %       corresponding to the underlying edges
            %
            %   v = qenergy(obj, s2, sst, mu);
            %
            %       This syntax returns v = (1/2) * tr(J * (C + mu*mu')).
            %
            
            jdv = obj.Jdv;
            w = obj.ew;
            
            if nargin >= 4
                s2 = s2 + mu .* mu;            
                sst = sst + mu(obj.es) .* mu(obj.et);                                
            end
            
            v = (jdv' * s2) / 2 + (w' * sst); 
        end
                        
        
        function v = qenergy_x(obj, x)
            % Compute the (1/2) * x' * J * x
            %
            %   v = qenergy_x(obj, x);
            %       computes the value (1/2) * x' * J * x. 
            %       In the input, x can be an nv x 1 vector or a matrix
            %       of nv x K. Then in the output, v is a vector of 
            %       size 1 x K, with v(k) corresponding to x(:,k).
            %
                        
            if ~(isfloat(x) && ndims(x) == 2 && size(x,1) == obj.nv)
                error('sgmrf:qenergy_x:invalidarg', ...
                    'x should be a numeric matrix with size(x,1) == nv');
            end
            
            if size(x,2) == 1
                xs = x(obj.es);
                xt = x(obj.et);
            else
                xs = x(obj.es, :);
                xt = x(obj.et, :);
            end            
            
            v = qenergy(obj, x.^2, xs .* xt);
        end

        
        function v = qenergy_sr(obj, sigma, rho)
            % Compute (1/2) * tr(J * C) using sigma and rho
            %
            %   v = qenergy_sigma_rho(obj, sigma, rho)
            %       computes the quadratic energy using marginal
            %       standard deviation given by sigma and the correlation
            %       coefficients given by rho
            %
            
            v = qenergy(obj, sigma .^ 2, ...
                rho .* (sigma(obj.es).*sigma(obj.et)) );
        end
        
    end
    
    
    
    methods(Static)
        
        function obj = amodel(g0, a)
            % Constructs an attractive gaussian mrf
            %
            %   obj = sgmrf.amodel(g0, a);
            %       constructs an attractive Gaussian MRF, of which 
            %       the information matrix is defined to be
            %
            %           J = laplacemat(g0, a);
            %
            %       Here, g0 should be an undirected gr_adjlist object,
            %       and a is either a scalar or an nv x 1 vector.
            %
            
            if ~(isa(g0, 'gr_adjlist') && g0.dtype == 'u')
                error('sgmrf:amodel:invalidarg', ...
                    'g0 should be an undirected gr_adjlist object.');
            end
            n = g0.nv;
            
            if ~( isfloat(a) && (isscalar(a) || isequal(size(a), [n 1])) )
                error('sgmrf:invalidarg', ...
                    'a should be either a scalar or an n x 1 numeric vector.');
            end
                        
            jdv = aggreg(g0.ew, n, g0.es+1, 'sum') + a;
            jgr = set_weights(g0, -g0.ew);                        
            obj = sgmrf(jdv, jgr);
        end
        
    end
    
    
end

