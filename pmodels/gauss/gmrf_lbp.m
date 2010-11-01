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
        
        
        function infer_Js(obj, vs)
            % Perform inference on information matrices
            %
            %   obj.infer_Js(vs);
            %       performs inference on information matrices for
            %       local models by updating messages sent from 
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

        
        function infer_hs(obj, vs)
            % Perform inference on potential vectors
            %
            %   obj.infer_hs(vs);
            %       performs inference on potential vectors for
            %       local models by updating messages sent from
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
    
    
end


