classdef gatrwa < handle
    % The class for tree-reweighted approximation of Gaussian model
    %
    
    % Created by Dahua Lin, on Feb 6, 2011
    %
    
    
    properties(GetAccess='public', SetAccess='private')
        nv;         % the number of nodes
        ne;         % the number of edges
        
        egraph;     % the underlying weighted edge graph for J
        Jdv;        % diagonal values of J
        eprob;      % the edge appearance probabilities [m x 1]
    end
       
    properties
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
            if ~isa(J, 'double'); J = double(J); end                
            
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
        
                
        function C = pseudo_cov(obj, s, r)
            % Get the pseudo covariance matrix from current solution
            %
            %   C = pseudo_cov(obj);
            %       returns the pseudo covariance matrix based on the
            %       current solution.
            %
            %   C = pseudo_cov(obj, s, r);
            %       returns the pseudo covariance matrix based on the
            %       solution given by s (sigma) and r (rho).
            %
            
            n = obj.nv;
            m = obj.ne;
            if nargin == 1
                s = obj.sigma;
                r = obj.rho;
            else
                if ~(isfloat(s) && isequal(size(s), [n 1]))
                    error('gatrwa:pseudo_cov:invalidarg', ...
                        's should be an n x 1 numeric vector.');
                end
                if ~(isfloat(r) && isequal(size(r), [m 1]))
                    error('gatrwa:pseudo_cov:invalidarg', ...
                        'r should be a 2m x 1 numeric vector.');
                end
            end
            
            eg = obj.egraph;            
            es = double(eg.es(1:m) + 1);
            et = double(eg.et(1:m) + 1);
            
            cvv = s .^ 2;
            cst = r .* (s(es) .* s(et));
            
            i = [(1:n)'; es; et];
            j = [(1:n)'; et; es];
            v = [cvv; cst; cst];
            
            C = sparse(i, j, v, n, n);            
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
        
        
        %% Inference steps
        
        function initialize(obj, s, r)
            % Initialize the solution
            %
            %   obj.initialize();
            %       This function initializes the solution of sigma 
            %       and rho.
            %           
            %       By default, sigma is initialized to be 1./sqrt(Jdv);
            %       and rho are set to all zeros.
            %
            %   obj.initialize(s, r);
            %       initializes the solution using s for sigma and 
            %       r for rho.
            %
            
            if nargin < 2
                obj.sigma = 1 ./ sqrt(obj.Jdv);
            else
                if ~(isfloat(s) && isequal(size(s), [obj.nv 1]))
                    error('gatrwa:initialize:invalidarg', ...
                        's should be an n x 1 numeric vector.');
                end
            end
            
            if nargin < 3
                obj.rho = zeros(obj.ne, 1);
            else
                if ~(isfloat(r) && isequal(size(r), [obj.ne 1]))
                    error('gatrwa:initialize:invalidarg', ...
                        'r should be a m x 1 numeric vector.');
                end
            end
            
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
            
            r = gatrwa_cimp(2, obj.egraph, obj.sigma, obj.eprob);
            obj.rho = r;           
        end
        
        
        function cupdate(obj, ord)
            % Perform combined updates
            %
            %   obj.cupdate(ord);
            %
            %       updates the sigma values of specified vertices along
            %       with the rho values of their incident edges
            %
            
            if ~(isvector(ord) && isnumeric(ord))
                error('gatrwa:invalidarg', 'ord should be a numeric vector.');
            end
            ord = int32(ord - 1);
                        
            [sig_, rho_] = gatrwa_cimp(3, obj.egraph, obj.Jdv, ...
                obj.sigma, obj.rho, obj.eprob, ord);
            
            obj.sigma = sig_;
            obj.rho = rho_;
        end

        
        %% Integrated solution
        
        function info = solve(obj, varargin)
            % solve through iterative cupdates
            %
            %   solve(obj, ...);
            %       solves the optimal solution through performing
            %       a series of combined updates until convergence.
            %
            %       One can specify further options to control the 
            %       process.
            %       - 'maxiter':    the maximum number of iterations
            %                       (default = 5000)
            %       - 'tolfun':     the tolerance of objective change
            %                       at convergence (default = 1e-8)
            %       - 'tolx':       the tolerance of sigma change at
            %                       convergence (default = 1e-6)
            %       - 'order':      the order of updating 
            %                       (default = 1 : nv)
            %       - 'display':    whether to display procedural info
            %                       (default = false)
            %       
            %   info = solve(obj, ...);
            %       it returns an info struct that contains the following
            %       fields:
            %       - 'niters':     the number of elapsed iterations
            %       - 'converged':  whether the procedure converged
            %       - 'energy':     the expected energy
            %       - 'entropy':    the approximated entropy
            %       - 'objv':       objective value: energy - entropy
            %                  
            
            % parse options
            
            maxiter = 5000;
            tolfun = 1e-8;
            tolx = 1e-6;
            order = [];
            display = false;
            
            if ~isempty(varargin)
                onames = varargin(1:2:end);
                ovals = varargin(2:2:end);
                
                if ~(length(onames) == length(ovals) && iscellstr(onames))
                    error('gatrwa:solve:invalidarg', ...
                        'The option list is invalid.');
                end
                
                for i = 1 : length(onames)
                    cn = onames{i};
                    cv = ovals{i};
                    switch lower(cn)
                        case 'maxiter'
                            if ~(isnumeric(cv) && isscalar(cv) && cv >= 1)
                                error('gatrwa:solve:invalidarg', ...
                                    'maxiter should be a positive integer.');
                            end
                            maxiter = cv;
                        case 'tolfun'
                            if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0)
                                error('gatrwa:solve:invalidarg', ...
                                    'tolfun should be a positive real scalar.');
                            end
                            tolfun = cv;
                        case 'tolx'
                            if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0)
                                error('gatrwa:solve:invalidarg', ...
                                    'tolx should be a positive real scalar.');
                            end
                            tolx = cv;
                        case 'order'
                            if ~(isnumeric(cv) && isvector(cv))
                                error('gatrwa:solve:invalidarg', ...
                                    'order should be a numeric vector.');
                            end
                            order = cv;
                        case 'display'
                            if ~(islogical(cv) && isscalar(cv))
                                error('gatrwa:solve:invalidarg', ...
                                    'display should be a logical scalar.');
                            end
                            display = cv;
                        otherwise
                            error('gatrwa:solve:invalidarg', ...
                                'Invalid option name %s', cn);
                    end
                end
            end
            
            if isempty(order)
                order = 1 : obj.nv;
            end
            
            % do solving
            
            it = 0;
            converged = false;
            
            energy = eval_energy(obj);
            entropy = eval_entropy(obj);
            objv = energy - entropy;
            
            if display
                fprintf('  Iter     obj.value     objv.ch    sigma.ch\n');
                fprintf('------------------------------------------------\n');
            end
            
            while ~converged && it < maxiter                
                it = it + 1;
                pre_objv = objv;
                pre_sig = obj.sigma;
                
                cupdate(obj, order);
                
                energy = eval_energy(obj);
                entropy = eval_entropy(obj);
                objv = energy - entropy;
                
                chf = objv - pre_objv;
                chx = max(abs(obj.sigma - pre_sig));
                converged = abs(chf) < tolfun && chx < tolx;
                
                if display
                    fprintf(' %5d  %12.4f  %10.3g  %10.3g\n', ...
                        it, objv, chf, chx);
                end
            end
            
            if display
                if converged
                    fprintf('converged with %d iterations.\n', it);
                else
                    fprintf('terminated without convergence.\n');
                end
            end            
            
            % output
            
            if nargout >= 1
                info = struct( ...
                    'niters', it, ...
                    'converged', converged, ...
                    'energy', energy, ...
                    'entropy', entropy, ...
                    'objv', objv);
            end            
        end        
        
        
        %% Objective function
        
        function [fv, gs, gr] = eval_energy(obj, s, r)
            % Evaluate the energy value
            %
            %   fv = obj.eval_energy();
            %   [fv, gs, gr] = obj.eval_energy();
            %
            %   fv = obj.eval_energy(s, r);
            %   [fv, gs, gr] = obj.eval_energy(s, r);
            %
            %       Evaluates the energy, which equals the expectation
            %       of (1/2) * (x-mu)' * J * (x-mu),
            %       under the pseudo-distribution given by the current
            %       solution.
            %
            %       If s and r are explicitly given, it uses the 
            %       solution given by s and r instead.
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
                    error('gatrwa:eval_energy:invalidarg', ...
                        's should be an n x 1 numeric vector.');
                end
                if ~(isfloat(r) && isequal(size(r), [m 1]))
                    error('gatrwa:eval_energy:invalidarg', ...
                        'r should be a m x 1 numeric vector.');
                end
            end
            
            eg = obj.egraph;
            
            jvv = obj.Jdv;
            jst = eg.ew(1:m);
                        
            lambda_v = s.^2;
            ss = s(eg.es + 1);
            sst = ss(1:m) .* ss(m+1:2*m);
            clear ss;
            
            lambda_e = r .* sst;  
                                                
            fv = (lambda_v' * jvv) / 2 + lambda_e' * jst;
            
            if nargout >= 3  % compute gradient                
                b = gatrwa_cimp(0, eg, s, r);                
                gs = jvv .* s + b;
                gr = jst .* sst;                
            end            
        end
        
        
        function [fv, gs, gr] = eval_entropy(obj, s, r)
            % Evaluate the energy value
            %
            %   fv = obj.eval_entropy();
            %   [fv, gs, gr] = obj.eval_energy();
            %
            %   fv = obj.eval_entropy(s, r);
            %   [fv, gs, gr] = obj.eval_energy(s, r);
            %
            %       Evaluates the approximated entropy, based on
            %       solution.
            %
            %       If s and r are explicitly given, it uses the 
            %       solution given by s and r instead.
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
                    error('gatrwa:eval_entropy:invalidarg', ...
                        's should be an n x 1 numeric vector.');
                end
                if ~(isfloat(r) && isequal(size(r), [m 1]))
                    error('gatrwa:eval_entropy:invalidarg', ...
                        'r should be a m x 1 numeric vector.');
                end
            end
            
            ep = obj.eprob;
            hc = (log(2*pi) + 1) / 2;            
            H = sum(log(s)) + hc * n;
            I = (-0.5) * (ep' * log(1 - r.^2));
            
            fv = H - I; 
            
            if nargout >= 3  % compute gradient
                gs = 1 ./ s;
                gr = - (ep .* (r ./ (1 - r.^2))); 
            end
            
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
            
            if nargin == 1
                s = obj.sigma;
                r = obj.rho;
            end
                        
            if nargout < 3                
                fv_e = eval_energy(obj, s, r);
                fv_h = eval_entropy(obj, s, r);
                fv = fv_e - fv_h;
            else
                [fv_e, gs_e, gr_e] = eval_energy(obj, s, r);
                [fv_h, gs_h, gr_h] = eval_entropy(obj, s, r);
                fv = fv_e - fv_h;
                gs = gs_e - gs_h;
                gr = gr_e - gr_h;
            end            
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
        
    end
    
    
end

