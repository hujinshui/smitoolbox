classdef gmrf_lbp < handle
    % The class for Loopy Belief propagation on Gaussian MRF
    %
    
    % History
    % -------
    %   - Created by Dahua Lin, on Nov 1, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')        
        nnodes;     % the number of nodes
        graph;      % the underlying undirected graph
        
        nb_nds;     % the cell array of neighbor nodes (n x 1)
        nb_oes;     % the cell array of outgoing neighbor edges (n x 1)
        nb_ies;     % the cell array of incoming neighbor edges (n x 1)
        
        dims;   % the dimension of each node
        Js;     % sub information matrices
        hs;     % sub potential vectors 
        Rs;     % (cross) information matrices (from parent to self)
        
        m_Ls;   % the L-matrices for generating messages
        m_Js;   % the messages for information matrices
        m_hs;   % the messages for potential vectors
        
        r_Js;   % the resultant information matrices of local models
        r_hs;   % the resultant potential vectors of local models        
    end
    
    
    methods
        
        function obj = gmrf_lbp(gm)
            % constructs a loopy belief-propagation solution on GMRF
            %
            %   obj = gmrf_lbp(gm);
            %       constructs a loopy belief-propagation solution on
            %       a Gaussian MRF gm, and make necessary initialization
            %       for the inference on information matrices.
            %
            
            if ~isa(gm, 'gaussmrf')
                error('gmrf_lbp:invalidarg', 'gm should be a gaussmrf object.');
            end
            
            % set graph and neighbor hood structure
            
            g = gm.graph;
            n = g.n;
            ne = g.m;
            
            nbs = cell(n, 1);
            oes = cell(n, 1);
            ies = cell(n, 1);
            
            for i = 1 : n
                d = g.o_ds(i);
                if d > 0
                    b = g.o_os(i);
                    v = g.o_ns(b + (1:d));
                    oe = g.o_es(b + (1:d));
                    ie = zeros(size(oe), 'int32');
                    ie(oe < ne) = oe(oe < ne) + ne;
                    ie(oe >= ne) = oe(oe >= ne) - ne;
                    
                    nbs{i} = v' + 1;
                    oes{i} = oe' + 1;
                    ies{i} = ie' + 1;
                end
            end
            
            obj.nnodes = n;
            obj.graph = g;  
            obj.nb_nds = nbs;
            obj.nb_oes = oes;
            obj.nb_ies = ies;
            
            % set basic info
            
            ds = gm.dims;            
            obj.dims = ds;                        
            obj.Js = gm.Js;
            obj.hs = [];
            obj.Rs = gm.Rs;
                
            obj.r_Js = cell(n, 1);
            obj.r_hs = [];        
            
            % initialize messages
            
            s = g.s + 1;
            t = g.t + 1;
            
            m = numel(s);
            mLs = cell(m, 1);
            mJs = cell(m, 1);            
            
            for k = 1 : m                
                di = ds(s(k));
                dj = ds(t(k));
                
                mLs{k} = zeros(dj, di);
                mJs{k} = zeros(dj, dj);                                
            end
            
            obj.m_Ls = mLs;
            obj.m_Js = mJs;
            obj.m_hs = [];
        end
        
        
        function r = infer_Js(obj, varargin)
            % perform inferences on information matrices
            %
            %   r = obj.infer_Js(...);
            %       performs inference on information matrices by 
            %       iteratively updating relevant messages.
            %
            %       One can specify additional options in form of
            %       name/value pairs to control the inference process:
            %       
            %       - 'order':  the order of nodes to be updated in 
            %                   each iteration. By default it is 1:n.
            %
            %       - 'maxiter':  the maximum number of iterations
            %                     (by default = 100)
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
            
            [order, maxiter, tol] = parse_infer_options(obj, varargin);
            n = obj.nnodes;
                        
            % initialize
            
            converged = false;
            it = 0;
            
            % iterations
            
            while ~converged && it < maxiter
                
                pre_Js = obj.r_Js;
                it = it + 1;                
                
                % do update
                update_Js(obj, order);
                
                % compare updates and determine convergence
                diffs = zeros(1, n);
                for i = 1 : n
                    pJ = pre_Js{i};
                    cJ = obj.r_Js{i};
                    
                    if ~isempty(pJ)
                        diffs(i) = Linfdiff(pJ, cJ);
                    else
                        diffs(i) = inf;
                    end                    
                end
                converged = max(diffs) < tol;                
                
            end
            
            % output
            if nargout >= 1
                r = struct('niters', it, 'converged', converged);
            end                        
        end
        
        
        
        function update_Js(obj, vs)
            % update messages about information matrices
            %
            %   obj.update_Js(vs);
            %       updates messages on information matrices sent from 
            %       a sequence of nodes specified by vs.
            %
            
            % verify input arguments
            
            if ~(isnumeric(vs) && isvector(vs))
                error('gmrf_lbp:infer_Js:invalidarg', ...
                    'vs should be a numeric vector.');
            end
            if size(vs, 1) > 1; vs = vs.'; end
            
            % prepare data
            
            nbs = obj.nb_nds;
            oes = obj.nb_oes;
            ies = obj.nb_ies;
            
            Js_ = obj.Js;
            Rs_ = obj.Rs;
            
            % update messages from each node
            for k = vs
                J = Js_{k};
                nnb = numel(nbs{k});
                
                if nnb > 0
                    % collect incoming messages
                    ie = ies{k};
                    if nnb == 1
                        J = J + obj.m_Js{ie};
                    else
                        J = J + sum(cat(3, obj.m_Js{ie}), 3);
                    end
                    
                    % set outgoing messages
                    oe = oes{k};
                    for j = 1 : nnb
                        ci = ie(j);
                        co = oe(j);
                        
                        J_uk = Rs_{ci};
                        L = J_uk / (J - obj.m_Js{ci});
                        
                        obj.m_Ls{co} = L;
                        obj.m_Js{co} = -L * J_uk';
                    end
                end
                
                % update marginal
                obj.r_Js{k} = J;
            end
            
        end
        
        
        function initialize_hs(obj, hs)
            % Initialize relevant structure & messages for inferring hs
            
            n = obj.nnodes;
            if ~(iscell(hs) && isvector(hs) && numel(hs) == n)
                error('gmrf_lbp:invalidarg', ...
                    'hs should be a cell vector with n cells.');
            end
            
            ds = obj.dims;
            K = size(hs{1}, 2);
            for i = 1 : n
                h = hs{i};
                if ~(isfloat(h) && isequal(size(h), [ds(i) K]))
                    error('gmrf_lbp:invalidarg', ...
                        'some h in hs is not valid.');
                end
            end
            
            obj.hs = hs;               
            obj.r_hs = cell(n, 1);
                        
            g = obj.graph;
            t = g.t + 1;
            m = numel(t);
            
            mhs = cell(m, 1);
            for i = 1 : m
                mhs{i} = zeros(ds(t(i)), K);
            end
            obj.m_hs = mhs;
            
        end
        
        
        function clear_hs(obj)
            % Clear the associated potential vectors and associated results
            
            obj.hs = [];
            obj.r_hs = [];
            obj.m_hs = [];
        end

        
        
        function r = infer_hs(obj, varargin)
            % perform inferences on potential vectors
            %
            %   r = obj.infer_hs(...);
            %       performs inference on potential vectors by 
            %       iteratively updating relevant messages.
            %
            %       One can specify additional options in form of
            %       name/value pairs to control the inference process:
            %       
            %       - 'order':  the order of nodes to be updated in 
            %                   each iteration. By default it is 1:n.
            %
            %       - 'maxiter':  the maximum number of iterations
            %                     (by default = 100)
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
            
            [order, maxiter, tol] = parse_infer_options(obj, varargin);
            n = obj.nnodes;
                        
            % initialize
            
            converged = false;
            it = 0;
            
            % iterations
            
            while ~converged && it < maxiter
                
                pre_hs = obj.r_hs;
                it = it + 1;                
                
                % do update
                update_hs(obj, order);
                
                % compare updates and determine convergence
                diffs = zeros(1, n);
                for i = 1 : n
                    ph = pre_hs{i};
                    ch = obj.r_hs{i};
                    
                    if ~isempty(ph)
                        diffs(i) = Linfdiff(ph, ch);
                    else
                        diffs(i) = inf;
                    end                    
                end
                converged = max(diffs) < tol;                
                
            end
            
            % output
            if nargout >= 1
                r = struct('niters', it, 'converged', converged);
            end      
            
        end
        
        
        
        function update_hs(obj, vs)
            % Perform inference on potential vectors
            %
            %   obj.infer_hs(vs);
            %       updates messages about potential vectors sent from
            %       a sequence of nodes specified by vs.
            %
            
            % verify input arguments
            
            if ~(isnumeric(vs) && isvector(vs))
                error('gmrf_lbp:infer_Js:invalidarg', ...
                    'vs should be a numeric vector.');
            end
            if size(vs, 1) > 1; vs = vs.'; end
            
            % prepare data
            
            nbs = obj.nb_nds;
            oes = obj.nb_oes;
            ies = obj.nb_ies;
            
            hs_ = obj.hs;
            
            % update messages from each node
            for k = vs
                h = hs_{k};
                nnb = numel(nbs{k});
                
                if nnb > 0
                    % collect incoming messages
                    ie = ies{k};
                    if nnb == 1
                        h = h + obj.m_hs{ie};
                    else
                        h = h + sum(cat(3, obj.m_hs{ie}), 3);
                    end
                    
                    % set outgoing messages
                    oe = oes{k};
                    for j = 1 : nnb
                        ci = ie(j);
                        co = oe(j);
                                                                        
                        L = obj.m_Ls{co};
                        obj.m_hs{co} = -L * (h - obj.m_hs{ci});
                    end
                end
                
                % update marginal
                obj.r_hs{k} = h;
            end
            
        end
        
    end
    
    
    methods(Access='private')
        
        function [order, maxiter, tol] = parse_infer_options(obj, options)
            % parse input options for inference 
            
            order = 1 : obj.nnodes;
            maxiter = 100;
            tol = 1e-6;
            
            if ~isempty(options)
                onames = options(1:2:end);
                ovals = options(2:2:end);
                
                nopts = numel(onames);
                if ~(numel(ovals) == nopts && iscellstr(onames))
                    error('gmrf_lbp:invalidopt', ...
                        'The option list is invalid.');
                end
                
                for i = 1 : nopts
                    ov = ovals{i};
                    switch onames{i}
                        case 'order'
                            if ~(isnumeric(ov) && ndims(ov)==2 && size(ov,1)==1)
                                error('gmrf_lbp:invalidopt', ...
                                    'order should be a numeric row vector.');
                            end
                            order = ov;
                        case 'maxiter'
                            if ~(isnumeric(ov) && isscalar(ov) && ov > 0)
                                error('gmrf_lbp:invalidopt', ...
                                    'maxiter should be a positive scalar.');
                            end
                            maxiter = ov;
                        case 'tol'
                            if ~(isfloat(tol) && isscalar(tol) && isreal(tol) && tol > 0)
                                error('gmrf_lbp:invalidopt', ...
                                    'tol should be a positive real scalar.');
                            end
                            tol = ov;
                    end
                end
            end            
        end
        
    end
            
    
end


