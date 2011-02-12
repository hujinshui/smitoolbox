classdef agmrf
    % The class to represent an attractive Gaussian MRF
    %
    %   Based on an undirected graph G = (V, E), the energy function 
    %   is defined to be
    %
    %         (1/2) * sum_{u in V} a_v (x_v - y_v)^2
    %       + (1/2) * sum_{{u,v} in E} w_{uv} (x_u - x_v)^2
    %
    
    % Created by Dahua Lin, on Feb 12, 2011
    %
    
    properties(GetAccess='public', SetAccess='private') 
        nv;         % the number of nodes
        ne;         % the number of (internal) edges
        
        g_int;      % the underlying graph   
        w_ext;      % the weights to y [nv x 1]                        
    end
    
    
    methods
       
        function obj = agmrf(g, a)
            % Construct an attractive GMRF
            %
            %   obj = agmrf(g, a);
            %       constructs an attractive Gaussian MRF as formalized
            %       above.
            %
            %       Inputs: 
            %       - g:    the underlying graph
            %       - a:    the external weight vector [nv x 1]
            %
            
            if ~(isa(g, 'gr_adjlist') && g.dtype == 'u')
                error('agmrf:invalidarg', ...
                    'g should be an undirected gr_adjlist object.');
            end         
            n = g.nv;
            
            if ~( isfloat(a) && (isscalar(a) || isequal(size(a), [n 1])) )
                error('agmrf:invalidarg', ...
                    'a should be either a scalar or an n x 1 numeric vector.');
            end
            if ~isa(a, 'double'); a = double(a); end
            
            obj.nv = n;
            obj.ne = g.ne;
            obj.g_int = g;
            obj.w_ext = a;                        
        end
        
        
        function gm = to_gmrf(obj)
            % Create a GMRF object from this model
            %
            %   gm = to_gmrf(obj);
            %
            
            g = obj.g_int;
            a = obj.w_ext;
            
            jdv = aggreg(g.ew, g.nv, g.es+1, 'sum') + a;
            jgr = set_weights(g, -g.ew);                        
            gm = gmrf(jdv, jgr);            
        end
        
       
        function J = information_matrix(obj)
            % Compute the information matrix
            %
            %   J = information_matrix(obj);
            %
            
            n = obj.nv;
            a = obj.w_ext;
            g = obj.g_int;
            
            jdv = aggreg(g.ew, n, g.es+1, 'sum') + a;
            i = [(1:n)'; double(g.es+1)];
            j = [(1:n)'; double(g.et+1)];
            v = [jdv; -g.ew];
            if ~isa(v, 'double')
                v = double(v);
            end
            
            J = sparse(i, j, v, n, n);
        end
        
    end
    
    
    methods
        
        %% Energy computation
        
        function v = ext_energy(obj, x, y)
            % Compute the external energy
            %
            %   v = obj.ext_energy(x, y);
            %
            
            a = obj.w_ext;
            e = (x - y).^2;
            v = ( a' * e ) / 2;            
        end
        
        
        function v = int_energy(obj, x)
            % Compute the internal energy
            %
            %   v = obj.int_energy(x, y);
            %
            
            g = obj.g_int;
            w = g.ew(1:g.ne);
            
            e = gmrf_cimp(2, g.ne, g.es, g.et, x);
            v = ( w' * e ) / 2;
        end
        
        function v = energy(obj, x, y)
            % Compute the energy given x and y
            %
            %   v = obj.energy(x, y);
            %       evaluates the energy value as defined above
            %
            
            v = ext_energy(obj, x, y) + int_energy(obj, x);
        end
        
        
        %% Expected energy computation
        
        function v = ext_energy_ep(obj, mx, vs, y)
            % Compute the expected external energy given distribution of x
            %
            %   v = obj.ext_energy_ep(obj, mx, vs, y);
            %       evaluates the expected external energy
            %       
            %   Inputs:           
            %       - mx:   the mean of x
            %       - vs:   the variance of x
            %       - y:    the external value 
            %
            
            a = obj.w_ext;
            e = (mx - y) .^ 2 + vs;
            v = ( a' * e ) / 2;            
        end
        
        
        function v = int_energy_ep(obj, mx, vs, cs)
            % Compute the expected internal energy given distribution of x
            %
            %   v = obj.int_energy_ep(obj, mx, vs, cs);
            %
            %   Inputs:
            %       - mx:   the mean of x
            %       - vs:   the variance of x
            %       - cs:   the covariance at edges
            %
            
            g = obj.g_int;
            w = g.ew(1:g.ne);
            
            e1 = gmrf_cimp(2, g.ne, g.es, g.et, mx);
            e2 = gmrf_cimp(0, g.ne, g.es, vs) + gmrf_cimp(0, g.ne, g.et, vs);           
            
            e = e1 + e2 - 2 * cs;
            v = ( w' * e ) / 2;            
        end
        
        
        function v = energy_ep(obj, mx, vs, cs, y)
            % Compute the expected energy given distribution of x
            %   v = obj.int_energy_ep(obj, mx, vs, cs);
            %
            %   Inputs:
            %       - mx:   the mean of x
            %       - vs:   the variance of x
            %       - cs:   the covariance at edges   
            %       - y:    the external values
            %
            
            v = ext_energy_ep(obj, mx, vs, y) + ...
                int_energy_ep(obj, mx, vs, cs);            
        end        
        
    end                
    
end

