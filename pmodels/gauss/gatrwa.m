classdef gatrwa < handle
    % The class for tree-reweighted approximation of Gaussian model
    %
    
    % Created by Dahua Lin, on Feb 6, 2011
    %
    
    properties
        
        nv;         % the number of nodes
        ne;         % the number of edges
        
        egraph;     % the underlying weighted edge graph for J
        Jdv;        % diagonal values of J
        eprob;      % the edge appearance probabilities [m x 1]
        
        sigma;      % the per-node sigma values [n x 1]: sqrt(var)
        rho;        % the per-edge coefficient correlation [m x 1]
    end
    
    
    methods
        
        function obj = gatrwa(J)
            % construct an object for tree-reweighted approximation
            %
            %   obj = gatrwa(J);
            %   obj = gatrwa(J, ep);
            %
            %       In the input, J is an information matrix. 
            %
            
            if ~(isfloat(J) && ndims(J) == 2 && size(J,1) == size(J,2))
                error('gatrwa:invalidarg', ...
                    'J should be a numeric symmetric matrix .');
            end
            
            % construct graph
            
            n = size(J, 1);
            [s, t, Jv] = find(J);
            
            si = find(s < t);
            s = s(si);
            t = t(si);
            Jv = Jv(si);
            
            eg = gr_adjlist.from_edges('u', n, s, t, Jv);
            
            obj.nv = n;
            obj.ne = eg.ne;
            obj.egraph = eg;
            obj.Jdv = full(diag(J));
            
            obj.eprob = [];
            obj.sigma = [];
            obj.rho = [];            
        end
        
        
        function set_eprob(obj, ep)
            % Set the edge appearance probability vector
            %         
            %   obj.set_eprob(ep);
            %
            %       Here, ep should be a vector of size m x 1.
            %       ep(i) is the value for the edge corresponding to 
            %       the i-th edge in obj.graph.
            %
            %       If one wants to set the same value to all edges,
            %       ep can be a scalar.
            %
            
            m = obj.ne;
            if ~(isfloat(ep) && (isscalar(ep) || isequal(size(ep), [m 1])))
                error('gatrwa:set_eprob:invalidarg', ...
                    'ep should be either a scalar or a numeric vector of size m x 1.');
            end
            
            if isscalar(ep);
                obj.eprob = constmat(m, 1, ep);
            else
                obj.eprob = ep;
            end
        end
        
        
        function initialize(obj)
            % Initialize the solution
            %
            %   obj.initialize();
            %       This function initializes the solution of sigma 
            %       and rho.
            %           
            %       By default, sigma is initialized to be 1./sqrt(Jdv);
            %       and rho are set to all zeros.
            
            obj.sigma = 1 ./ sqrt(obj.Jdv);
            obj.rho = zeros(obj.ne, 1);
        end
        
        
        function update_sigma(obj, ord)
            % Update the sigma values
            %
            %   obj.update_sigma();
            %   obj.update_sigma(ord);
            %
            %       Updates the sigma values in specified order. 
            %       If ord is omitted, it is equivalent to set ord to 1:n.
            %  
            
            if nargin < 2
                ord = int32(0:obj.nv-1);
            else
                if ~(isvector(ord) && isnumeric(ord))
                    error('gatrwa:invalidarg', 'ord should be a numeric vector.');
                end
                ord = int32(ord - 1);
            end
                                                
            sig_ = gatrwa_cimp(1, obj.egraph, obj.Jdv, obj.sigma, obj.rho, ord);
            obj.sigma = sig_;
        end        
        
        
        function update_rho(obj)
            % Update all rho values
            %
            %   obj.update_rho();
            %
            %       Updates all rho values. Since the updating of each
            %       rho value is independent, the user do not need to
            %       set the order as in the updating of sigma.
            %
            
            m = obj.ne;
            eg = obj.egraph;                                         
            ep = obj.eprob;
            
            % compute a
            ss = obj.sigma(eg.es + 1);
            a = eg.ew(1:m) .* (ss(1:m) .* ss(m+1:2*m));
            
            r = (ep - sqrt(ep.^2 + 4 * a.^2)) ./ (2 * a);
            obj.rho = r;           
        end
        
        
        function f = objfunc(obj)
            % Get the objective function
            %
            %   f = obj.objfunc();
            %
            
            function [v, g] = gatrwaf(x)
                
                n = obj.nv;
                m = obj.ne;
                
                s = x(1:n);
                r = x(n+1:n+m);
                
                if nargout <= 1
                    v = eval_objv(obj, s, r);
                else
                    [v, gs, gr] = eval_objv(obj, s, r);
                    g = [gs; gr];
                end                
            end
            
            f = @gatrwaf;
        end
        
        
        
        function [fv, gs, gr] = eval_objv(obj, s, r)
            % Evaluate the objective value of a solution
            %
            %   fv = obj.eval_objv();
            %   [fv, gs, gr] = obj.eval_objv();
            %
            %   fv = obj.eval_objv(s, r);
            %   [fv, gs, gr] = obj.eval_objv(s, r);
            %
            %       evaluate the objective value for the given solution.
            %       Here, s is an n x 1 vector for sigma, and r is 
            %       an m x 1 vector for rho.
            %
            %       If s and r are omitted, the current solution would 
            %       be used.
            %
            %       The function returns the objective value fv, and if
            %       the second and third output are required, it also 
            %       returns the corresponding derivatives, via gs and gr.            
            %       
            
            n = obj.nv;
            m = obj.ne;
            
            if nargin == 1
                s = obj.sigma;
                r = obj.rho;
            else
                if ~(isfloat(s) && isequal(size(s), [n 1]))
                    error('gatrwa:eval_objv:invalidarg', ...
                        's should be an n x 1 numeric vector.');
                end
                if ~(isfloat(r) && isequal(size(r), [m 1]))
                    error('gatrwa:eval_objv:invalidarg', ...
                        'r should be a 2m x 1 numeric vector.');
                end
            end
            
            eg = obj.egraph;
            ep = obj.eprob;
            jvv = obj.Jdv;
            jst = eg.ew(1:m);
            hc = (log(2*pi) + 1) / 2;
            
            lambda_v = s.^2;
            ss = s(eg.es + 1);
            sst = ss(1:m) .* ss(m+1:2*m);
            clear ss;
            
            lambda_e = r .* sst;  
                                                
            fv0 = (lambda_v' * jvv) / 2 + lambda_e' * jst;
            H = sum(log(s)) + hc * n;
            I = (-0.5) * (ep' * log(1 - r.^2));
            
            fv = fv0 - H + I;
            
            if nargout >= 3  % compute gradient
                
                b = gatrwa_cimp(0, eg, s, r);
                
                gs = jvv .* s + b - 1 ./ s;
                gr = jst .* sst + ep .* (r ./ (1 - r.^2));                
            end
            
        end        
    end
    
    
end

