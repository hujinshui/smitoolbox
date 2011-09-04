classdef gmrf_blk_tbp < handle
    % The class for belief-propagation on block-tree-structured Gaussian MRF
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
        
        function obj = gmrf_blk_tbp(gm, varargin)
            % constructs a tree belief-propagation solution
            %
            %   obj = gmrf_blk_tbp(gm, vs, ps, es);
            %       constructs a tree belief-propagation solution for
            %       the Gaussian MRF gm. 
            %
            %       The computation tree is specified by vs, ps, and es. 
            %       Here, vs is the topological order (downward order
            %       from root), ps are the parents corresponding to vs,
            %       es are the corresponding downward edge indices.
            %
            %   obj = gmrf_blk_tbp(gm, seeds);
            %       uses the root seeds for constructing the computation 
            %       tree.
            %       
            
            % verify input
            
            if ~isa(gm, 'gmrf_blk')
                error('gmrf_blk_tbp:invalidarg', ...
                    'gm should be a gmrf_blk object.');
            end
                        
            if nargin == 2
                seeds = varargin{1};
                [vs, ps, es, tf] = gr_bfs_trees(gm.graph, seeds);
                if ~tf
                    error('gmrf_blk_tbp:invalidarg', ...
                        'The graph of gm is not in tree-structure.');
                end
            else            
                vs = varargin{1};
                ps = varargin{2};
                es = varargin{3};
                
                if ~(isnumeric(vs) && isnumeric(ps) && isnumeric(es) ...
                        && isvector(vs) && isequal(size(vs), size(ps), size(es)))
                    error('gmrf_blk_tbp:invalidarg', ...
                        'vs, ps, and es should be numeric vectors of the same size.');
                end
            end
            
            % construct the tree
            
            n = gm.nnodes;
            parents = zeros(n, 1);
            parents(vs) = ps;
            T = gr_tree.from_parents(parents);
            
            % prepare data structure for results
            
            Js_ = gm.Js;
            m = numel(es);
            gRs = gm.Rs;
            Rs_ = cell(m, 1);
            for k = 1 : m
                e = es(k);
                if e > 0
                    Rs_{vs(k)} = gRs{e};
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
            
            n = obj.tree.nv;
            if ~(iscell(hs) && isvector(hs) && numel(hs) == n)
                error('gmrf_blk_tbp:invalidarg', ...
                    'hs should be a cell vector with n cells.');
            end
            
            ds = obj.dims;
            K = size(hs{1}, 2);
            for i = 1 : n
                h = hs{i};
                if ~(isfloat(h) && isequal(size(h), [ds(i) K]))
                    error('gmrf_blk_tbp:invalidarg', ...
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
                v = vs(i);
                
                % collect messages from the children of v to v 
                J = Js_{v};
                if ~T.is_leaf(v)
                    cs = T.children_of(v);
                    if numel(cs) == 1
                        J = J + obj.up_Js{cs};
                    else
                        J = J + sum(cat(3, obj.up_Js{cs}), 3);
                    end
                end
                obj.r_Js{v} = J;
                
                % set messages to parent (if not root)
                if ~T.is_root(v)
                    Jpv = obj.Rs{v};
                    L = Jpv / J;
                    uJ = -(L * Jpv');
                    
                    obj.up_Ls{v} = L;
                    obj.up_Js{v} = uJ;
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
                v = vs(i);
                
                % collect messages from the children of v to v
                h = hs_{v};
                if ~T.is_leaf(v)
                    cs = T.children_of(v);
                    if numel(cs) == 1
                        h = h + obj.up_hs{cs};
                    else
                        h = h + sum(cat(3, obj.up_hs{cs}), 3);
                    end
                end
                obj.r_hs{v} = h;
                
                % set messages to parent (if not root)
                if ~T.is_root(v)
                    L = uLs{v};
                    obj.up_hs{v} = -L * h;
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
                if ~T.is_root(p)
                    J = J + obj.dn_Js{p};
                    obj.r_Js{p} = J;
                end
                
                % set messages to children (if not leaf)
                if ~T.is_leaf(p)
                    cs = T.children_of(p);
                    cs = cs(:).';
                    for c = cs                        
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
                if ~T.is_root(p)
                    h = h + obj.dn_hs{p};
                    obj.r_hs{p} = h;
                end
                
                % set messages to children (if not leaf)
                if ~T.is_leaf(p)
                    cs = T.children_of(p);
                    cs = cs(:).';
                    for c = cs                       
                        L = dLs{c};
                        
                        dh = -L * (h - uhs{c});
                        obj.dn_hs{c} = dh;
                    end
                end
            end            
        end
                               
    end
        
end


