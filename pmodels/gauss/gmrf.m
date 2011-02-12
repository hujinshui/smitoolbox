classdef gmrf
    % The class to represent a Gaussian MRF 
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
    end
        
    methods
        
        function obj = gmrf(v1, v2)
            % Construct a Gaussian MRF
            %
            %   obj = gmrf(J);
            %       constructs a Gaussian MRF from the information matrix.
            %       
            %   obj = gmrf(jdv, jg);
            %       constructs a Gaussian MRF from the diagonal values
            %       in jdv, and an off-diagonal graph jg.
            %   
            
            if nargin == 1
                J = v1;                
                if ~(isfloat(J) && ndims(J) == 2 && size(J,1) == size(J,2))
                    error('gmrf:invalidarg', 'J should be a symmetric matrix.');
                end
                
                jdv = full(diag(J));
                [s, t, w] = find(J);
                i = find(s < t);
                jgr = gr_adjlist.from_edges('u', size(J,1), s(i), t(i), w(i));
                
            else
                jdv = v1;
                jgr = v2;
                
                if ~(isa(jgr, 'gr_adjlist') && jgr.dtype == 'u')
                    error('gmrf:invalidarg', 'jgr should be an undirected graph.');
                end
                
                if ~(isfloat(jdv) && isequal(size(jdv), [jgr.nv, 1]))
                    error('gmrf:invalidarg', 'jdv should be a vector of size nv x 1.');
                end
            end
            
            obj.nv = jgr.nv;
            obj.ne = jgr.ne;
            
            obj.Jdv = jdv;
            obj.Jgr = jgr;    
        end
        
        
        function J = information_matrix(obj)
            % Get the information matrix J
            %
            %   J = information_matrix(obj);
            %
            % Note the returned matrix J is an nv x nv sparse matrix.
            %
                                    
            n = obj.nv;                 
            g = obj.Jgr;
            
            i = [(1:n)'; double(g.es+1)];
            j = [(1:n)'; double(g.et+1)];
            v = [obj.Jdv; g.ew];
            if ~isa(v, 'double')
                v = double(v);
            end
            
            J = sparse(i, j, v, n, n);
        end
        
                             
        function v = energy(obj, v, c)
            % Compute the energy value tr(J * C) / 2
            %
            %   v = qenergy(obj, v, c);
            %
            %       It returns the value as computed below
            %
            %       v = <Jdv, v> / 2 + <Jgr.ew, c>
            %
            %       This value equals tr(J * C) / 2, if C is a covariance
            %       matrix such that v captures its diagonal values, and
            %       c gives the covariance at graph edges.
            %
                                 
            v = (obj.Jdv' * v) / 2 + (obj.Jgr.ew(1:obj.ne)' * c); 
        end
                   
        
        function [vs, cs] = sr_to_vcs(obj, sigma, rho)
            % Converts (sigma, rho) to (vs, cs)
            %
            %   [vs, cs] = sr_to_vcs(obj, sigma, rho);
            %
            
            g = obj.Jgr;
            vs = sigma.^2;
            cs = rho .* gmrf_cimp(1, g.ne, g.es, g.et, sigma);            
        end
        
        
        
        function v = energy_sr(obj, sigma, rho)
            % Compute (1/2) * tr(J * C) using sigma and rho
            %
            %   v = qenergy_sigma_rho(obj, sigma, rho)
            %       computes the quadratic energy using marginal
            %       standard deviation given by sigma and the correlation
            %       coefficients given by rho
            %
            
            [vs, cs] = sr_to_vcs(obj, sigma, rho);            
            v = energy(obj, vs, cs);
        end
        
        
        function v = qterm(obj, x)
            % Compute the quadratic term (1/2) * x' * J * x
            %
            %   v = qterm(obj, x);
            %
            %       computes the value (1/2) * x' * J * x. 
            %       In the input, x can be an nv x 1 vector or a matrix
            %       of nv x K. Then in the output, v is a vector of 
            %       size 1 x K, with v(k) corresponding to x(:,k).
            %
                        
            if ~(isfloat(x) && ndims(x) == 2 && size(x,1) == obj.nv)
                error('gmrf:qenergy_x:invalidarg', ...
                    'x should be a numeric matrix with size(x,1) == nv');
            end
            if ~isa(x, 'double'); x = double(x); end
            
            g = obj.Jgr;
            xst = gmrf_cimp(1, g.ne, g.es, g.et, x);
            
            v = energy(obj, x.^2, xst);
        end

    end
    
    
    
    methods(Static)
        
        function obj = amodel(g0, a)
            % Constructs an attractive gaussian mrf
            %
            %   obj = gmrf.amodel(g0, a);
            %       constructs an attractive Gaussian MRF, of which 
            %       the information matrix is defined to be
            %
            %           J = laplacemat(g0, a);
            %
            %       Here, g0 should be an undirected gr_adjlist object,
            %       and a is either a scalar or an nv x 1 vector.
            %
            
            if ~(isa(g0, 'gr_adjlist') && g0.dtype == 'u')
                error('gmrf:amodel:invalidarg', ...
                    'g0 should be an undirected gr_adjlist object.');
            end
            n = g0.nv;
            
            if ~( isfloat(a) && (isscalar(a) || isequal(size(a), [n 1])) )
                error('gmrf:invalidarg', ...
                    'a should be either a scalar or an n x 1 numeric vector.');
            end
                        
            jdv = aggreg(g0.ew, n, g0.es+1, 'sum') + a;
            jgr = set_weights(g0, -g0.ew);                        
            obj = gmrf(jdv, jgr);
        end
        
    end
    
    
end

