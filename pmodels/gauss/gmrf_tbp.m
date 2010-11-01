classdef gmrf_tbp < handle
    % The class for belief-propagation on tree-structured Gaussian MRF
    %
    
    % Created by Dahua Lin, on Oct 31, 2010
    %
    
    
    properties(GetAccess='public', SetAccess='private')
        
        tree;   % the computation tree
        tord;   % the topological order
        
        dims;   % the dimension of each node
        Js;     % sub information matrices
        hs;     % sub potential vectors 
        Rs;     % (cross) information matrices (from parent to self)
        
        up_Ls;
        up_Js;  % the upward messages for information matrices (attached to children)
        up_hs;  % the upward messages for potential vectors (attached to children)
        
        dn_Ls;
        dn_Js;  % the downward messages for information matrices (attached to children)
        dn_hs;  % the downward messages for potential vectors (attached to children)
        
        r_Js;   % resultant information matrices
        r_hs;   % resultant potential vectors
    end
    
    
    methods
        
        function obj = gmrf_tbp(gm, varargin)
            % constructs a tree belief-propagation solution
            %
            %   obj = gmrf_tbp(gm, vs, ps, es);
            %       constructs a tree belief-propagation solution for
            %       the Gaussian MRF gm. 
            %
            %       The computation tree is specified by vs, ps, and es. 
            %       Here, vs is the topological order (downward order
            %       from root), ps are the parents corresponding to vs,
            %       es are the corresponding downward edge indices.
            %
            %   obj = gmrf_tbp(gm, seeds);
            %       uses the root seeds for constructing the computation 
            %       tree.
            %       
            
            % verify input
            
            if ~isa(gm, 'gaussmrf')
                error('gmrf_tbp:invalidarg', 'gm should be a gaussmrf object.');
            end
                        
            if nargin == 2
                seeds = varargin{1};
                [vs, ps, es, tf] = gr_bfs_trees(gm.graph, seeds);
                if ~tf
                    error('gmrf_tbp:invalidarg', ...
                        'The graph of gm is not in tree-structure.');
                end
            else            
                vs = varargin{1};
                ps = varargin{2};
                es = varargin{3};
                
                if ~(isnumeric(vs) && isnumeric(ps) && isnumeric(es) ...
                        && isvector(vs) && isequal(size(vs), size(ps), size(es)))
                    error('gmrf_bp:invalidarg', ...
                        'vs, ps, and es should be numeric vectors of the same size.');
                end
            end
            
            % construct the tree
            
            n = gm.nnodes;
            parents = zeros(n, 1);
            parents(vs) = ps;
            T = gr_tree(parents);
            
            % prepare data structure for results
            
            Js_ = gm.Js;
            m = numel(es);
            gRs = gm.Rs;
            Rs_ = cell(m, 1);
            for k = 1 : m
                e = es(k);
                if e > 0
                    Rs_{k} = gRs{e};
                end
            end
                        
            % set fields
            
            obj.tree = T;
            obj.tord = vs;
            
            obj.dims = gm.dims;
            obj.Js = Js_;
            obj.hs = [];
            obj.Rs = Rs_;
            
            obj.up_Ls = cell(n, 1);
            obj.up_Js = cell(n, 1);
            obj.up_hs = [];
            obj.dn_Ls = cell(n, 1);
            obj.dn_Js = cell(n, 1);
            obj.dn_hs = [];  
            
            obj.r_Js = cell(n, 1);
            obj.r_hs = [];            
        end
        
        
        function infer_Js(obj)
            % Inferring information matrices of all local models
            
            bottom_up_Js(obj);
            top_down_Js(obj);
        end          
        
        
        function initialize_hs(obj, hs)
            % Initialize relevant data structure for inferring hs
            
            n = obj.tree.n;
            if ~(iscell(hs) && isvector(hs) && numel(hs) == n)
                error('gmrf_tbp:invalidarg', ...
                    'hs should be a cell vector with n cells.');
            end
            
            ds = obj.dims;
            K = size(hs{1}, 2);
            for i = 1 : n
                h = hs{i};
                if ~(isfloat(h) && isequal(size(h), [ds(i) K]))
                    error('gmrf_tbp:invalidarg', ...
                        'some h in hs is not valid.');
                end
            end
            
            obj.hs = hs;            
        end
        
        
        function infer_hs(obj)
            % Inferring potential vectors of all local models
            
            bottom_up_hs(obj);
            top_down_hs(obj);
        end
        
        
        function clear_hs(obj)
            % Clear the associated potential vectors and associated results
            
            obj.hs = [];
            obj.r_hs = [];
            obj.up_hs = [];
            obj.dn_hs = [];            
        end
                            
    end
    
    
    
    methods        
        function bottom_up_Js(obj)
            % Performs Bottom-up collection for inferring Js along the tree
            
            T = obj.tree;
            Js_ = obj.Js;   
            
            vs = obj.tord;
            n = length(vs);
            
            % for each node in a bottom-up order (rev topological order)
            for i = n : -1 : 1
                c = vs(i);
                
                % collect messages from the children of c to c
                nc = T.ncs(c);                
                J = Js_{c};
                if nc > 0
                    cs = T.cs(T.os(c)+(1:nc)) + 1;
                    if nc == 1
                        J = J + obj.up_Js{cs};
                    else
                        J = J + sum(cat(3, obj.up_Js{cs}), 3);
                    end
                end
                obj.r_Js{c} = J;
                
                % set messages to parent (if not root)
                p = T.ps(c) + 1;
                if p > 0
                    Jpc = obj.Rs{c};
                    L = Jpc / J;
                    uJ = -(L * Jpc');
                    
                    obj.up_Ls{c} = L;
                    obj.up_Js{c} = uJ;
                end                
            end                        
        end
        
        
        function bottom_up_hs(obj)
            % Performs bottom-up collection for inferring hs along tree
            
            T = obj.tree;
            hs_ = obj.hs;
            uLs = obj.up_Ls;
            
            vs = obj.tord;
            n = length(vs);
            
            % for each node in a bottom-up order (rev topological order)
            for i = n : -1 : 1
                c = vs(i);
                
                % collect messages from the children of c to c
                h = hs_{c};
                nc = T.ncs(c);
                if nc > 0
                    cs = T.cs(T.os(c) + (1:nc)) + 1;
                    if nc == 1
                        h = h + obj.up_hs{cs};
                    else
                        h = h + sum(cat(3, obj.up_hs{cs}), 3);
                    end
                end
                obj.r_hs{c} = h;
                
                % set messages to parent (if not root)
                p = T.ps(c) + 1;
                if p > 0
                    L = uLs{c};
                    obj.up_hs{c} = -L * h;
                end
            end            
        end
        
        
        
        function top_down_Js(obj)
            % Performs top-down dispatch for inferring Js along the tree
            
            T = obj.tree;
            Rs_ = obj.Rs;
            uJs = obj.up_Js;
            
            vs = obj.tord;
            n = length(vs);
            
            % for each node in a top-down order (topological order)
            for i = 1 : n
                p = vs(i);
                
                % incorporate top-down message from parent of p to p
                J = obj.r_Js{p};
                pp = T.ps(p) + 1;
                if pp > 0
                    J = J + obj.dn_Js{p};
                    obj.r_Js{p} = J;
                end
                
                % set messages to children (if not leaf)
                nc = T.ncs(p);
                if nc > 0
                    cs = T.cs(T.os(p)+(1:nc)) + 1;
                    for j = 1 : nc
                        c = cs(j);
                        Jpc = Rs_{c};
                        L = ((J - uJs{c}) \ Jpc)';
                        dJ = -L * Jpc;
                        
                        obj.dn_Ls{c} = L;
                        obj.dn_Js{c} = dJ;
                    end
                end
            end            
        end
        
        
        function top_down_hs(obj)
            % Performs top-down dispatch for inferring the hs along tree
            
            T = obj.tree;
            dLs = obj.dn_Ls;
            uhs = obj.up_hs;
            
            vs = obj.tord;
            n = length(vs);
            
            % for each node in a top-down order (topological order)
            for i = 1 : n
                p = vs(i);
                
                % incorporate top-down message from parent of p to p
                h = obj.r_hs{p};
                pp = T.ps(p) + 1;
                if pp > 0
                    h = h + obj.dn_hs{p};
                    obj.r_hs{p} = h;
                end
                
                % set messages to children (if not leaf)
                nc = T.ncs(p);
                if nc > 0
                    cs = T.cs(T.os(p)+(1:nc)) + 1;
                    for j = 1 : nc
                        c = cs(j);
                        L = dLs{c};
                        
                        dh = -L * (h - uhs{c});
                        obj.dn_hs{c} = dh;
                    end
                end
            end            
        end
                               
    end
        
end


