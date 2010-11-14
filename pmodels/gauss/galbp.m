classdef galbp < handle
    % The class for loopy belief-propagation of Gaussian distribution
    %
    %   The class implements a scalar-level Gaussian LBP algorithm.
    %
    
    % Created by Dahua Lin, on Nov 3, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        J;          % the information matrix
        nnodes;     % the number of nodes (n)
        nedges;     % the number of edges (m)
        
        es;         % the source of edges (2m x 1)
        et;         % the target of edges (2m x 1)
        ew;         % the edge weights (2m x 1)
        rev_e;      % the corresponding indices of reversed edge (2m x 1)
        
        Jdv;        % the diagonal entries of J
        
        m_L;        % the L values of messages (2m x 1)
        m_J;        % the J messages (2m x 1)
        
        r_Js;       % the resultant vector of local precision (n x 1)                        
    end
    
    
    methods
        
        function obj = galbp(J)
            % constructs a Gaussian LBP solution 
            %
            %   obj = galbp(J);
            %       constructs a loopy BP solution on a Gaussian model
            %       whose information matrix is J.
            %
            
            % verify input arguments
            n = size(J, 1);
            if ~(isfloat(J) && ndims(J) == 2 && n == size(J,2) )
                error('galbp:invalidarg', 'J should be a numeric square matrix.');
            end
            
            % extract graph edges
                                                
            Jdv_ = full(diag(J));
            [s, t, w] = find(J);
            se = s < t;
            s = s(se);
            t = t(se);
            w = w(se);
            m = numel(s);                        
                        
            % set fields                        
            
            obj.nnodes = n;
            obj.nedges = m;
            
            obj.es = [s; t];
            obj.et = [t; s];
            obj.ew = [w; w];
            obj.rev_e = [m+1:2*m, 1:m].';
            obj.Jdv = Jdv_;
            
            obj.m_L = [];
            obj.m_J = zeros(2*m, 1);
            
            obj.r_Js = [];
            
        end
        
        
        function r = infer_Js(obj, varargin)
            % Perform inference of local information matrix
            %
            %   r = obj.infer_Js( ... );
            %       infers the local precision value (local information
            %       matrix).
            %
            %       One can specify additional options in form of
            %       name/value pairs to control the inference process:
            %       
            %       - 'maxiter':  the maximum number of iterations
            %                   (by default = 100)
            %
            %       - 'tol':    the tolerance of changes at convergence
            %                   the change is measured by the maximum 
            %                   L-inf norm of local information matrix 
            %                   difference. (by default = 1e-6)
            %
            %       The output is a struct recording the runtime
            %       information, which comprises two fields:
            %
            %       - 'niters':     the number of elapsed iterations
            %       - 'converged':  whether the iteration converges
            %
            
            % verify and parse input 
            
            [maxiter, tol] = parse_infer_options(obj, varargin);
            n = obj.nnodes;
                                                            
            % initialize                                    
            
            s = obj.es;
            t = obj.et;
            w = obj.ew;
            re = obj.rev_e;
            
            Jd = obj.Jdv;
            mJ = obj.m_J;
            
            converged = false;
            it = 0;
            
            % iterations
            
            while ~converged && it < maxiter
                
                pre_Jr = obj.r_Js;
                it = it + 1;                
                
                % do parallel update
                
                Jr = Jd + aggreg(mJ, n, t, 'sum');
                L = w ./ (Jr(s) - mJ(re));
                mJ = -L .* w;                
                
                % compare updates and determine convergence
                if ~isempty(pre_Jr)
                    converged = max(abs(Jr - pre_Jr)) < tol;
                end
            end
            
            % store the results
            
            obj.m_L = L;
            obj.m_J = mJ;
            obj.r_Js = Jr;
            
            
            % output
            if nargout >= 1
                r = struct('niters', it, 'converged', converged);
            end            
            
        end
        
        
        
        function [hr, r] = infer_hs(obj, h, varargin)
            % Perform inference of local potential vector
            %
            %   [hr, r] = obj.infer_hs(h, ... );
            %       infers the local potentials based on potential 
            %       vector h.
            %
            %       In output, hr is the inferred vector of local
            %       potentials.
            %
            %       One can specify additional options in form of
            %       name/value pairs to control the inference process:
            %       
            %       - 'maxiter':  the maximum number of iterations
            %                   (by default = 100)
            %
            %       - 'tol':    the tolerance of changes at convergence
            %                   the change is measured by the maximum 
            %                   L-inf norm of local information matrix 
            %                   difference. (by default = 1e-6)
            %
            %       The 2nd output is a struct recording the runtime
            %       information, which comprises two fields:
            %
            %       - 'niters':     the number of elapsed iterations
            %       - 'converged':  whether the iteration converges
            %
            
            % verify and parse input 
            
            [maxiter, tol] = parse_infer_options(obj, varargin);
            n = obj.nnodes;
            m = obj.nedges;
                                                            
            % initialize                                    
            
            s = obj.es;
            t = obj.et;
            re = obj.rev_e;
            
            L = obj.m_L;
            mh = zeros(2*m, 1);
            
            converged = false;
            it = 0;
            
            % iterations
            
            while ~converged && it < maxiter
                
                it = it + 1;                
                
                % do parallel update
                
                hr = h + aggreg(mh, n, t, 'sum');
                mh = -L .* (hr(s) - mh(re)); 
                
                % compare updates and determine convergence
                if it > 1
                    converged = max(abs(hr - pre_hr)) < tol;
                end
                pre_hr = hr;
                
            end
                        
            
            % output
            if nargout >= 1
                r = struct('niters', it, 'converged', converged);
            end            
            
        end        
        
    end
    
    
    
    methods(Access='private')
        
        function [maxiter, tol] = parse_infer_options(obj, options) %#ok<MANU>
            % parse input options for inference 
            
            maxiter = 100;
            tol = 1e-8;
            
            if ~isempty(options)
                onames = options(1:2:end);
                ovals = options(2:2:end);
                
                nopts = numel(onames);
                if ~(numel(ovals) == nopts && iscellstr(onames))
                    error('galbp:invalidopt', ...
                        'The option list is invalid.');
                end
                
                for i = 1 : nopts
                    ov = ovals{i};
                    switch onames{i}
                        case 'maxiter'
                            if ~(isnumeric(ov) && isscalar(ov) && ov > 0)
                                error('galbp:invalidopt', ...
                                    'maxiter should be a positive scalar.');
                            end
                            maxiter = ov;
                        case 'tol'
                            if ~(isfloat(tol) && isscalar(tol) && isreal(tol) && tol > 0)
                                error('galbp:invalidopt', ...
                                    'tol should be a positive real scalar.');
                            end
                            tol = ov;
                    end
                end
            end            
        end
        
    end
    
    
end


