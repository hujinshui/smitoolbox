classdef gmrf_tbp < handle
    % The class for belief-propagation on tree-structured Gaussian MRF
    %
    
    % Created by Dahua Lin, on Oct 31, 2010
    %
    
    
    properties(GetAccess='public', SetAccess='private')
        
        tree;   % the computation tree
        tord;   % the topological order
        
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
        
        
        function bottom_up_J(obj)
            % Performs Bottom-up collection along the tree
            
            vs = obj.tord;
            n = length(vs);
            
            for i = n : -1 : 1
                up_collect_J(obj, vs(i));
            end            
            
        end
        
        
        function top_down_J(obj)
            % Performs top-down dispatch along the tree
            
            vs = obj.tord;
            n = length(vs);
            
            for i = 1 : n
                down_dispatch_J(obj, vs(i));               
            end            
        end
        
        
        function up_collect_J(obj, c) 
            % compute upward message from a child node
            
            T = obj.tree;
            nc = T.ncs(c);
            
            J = obj.Js{c};
            if nc > 0
                cs = T.cs(T.os(c)+(1:nc)) + 1;
                if nc == 1
                    J = J + obj.up_Js{cs};
                else
                    J = J + sum(cat(obj.up_Js{cs}, 3), 3);
                end
            end
            obj.r_Js{c} = J;
            
            p = T.ps(c) + 1;
            if p > 0            
                Jpc = obj.Rs{c};
                L = Jpc / J;            
                uJ = -(L * Jpc');
                        
                obj.up_Ls{c} = L;
                obj.up_Js{c} = uJ;
            end
        end
        
        
        function down_dispatch_J(obj, p)
            % compute downward message from a parent node
            
            T = obj.tree;
            J = obj.r_Js{p};            
            Rs_ = obj.Rs;
            uJs = obj.up_Js;
            
            pp = T.ps(p) + 1;                        
            if pp > 0
                J = J + obj.dn_Js{p};
                obj.r_Js{p} = J;
            end
            
            nc = T.ncs(p);
            if nc > 0                
                cs = T.cs(T.os(p)+(1:nc)) + 1;
                for i = 1 : nc
                    c = cs(i);
                    Jpc = Rs_{c};
                    L = ((J - uJs{c}) \ Jpc)';
                    dJ = -L * Jpc;
                    
                    obj.dn_Ls{c} = L;
                    obj.dn_Js{c} = dJ;                    
                end                
            end            
        end
        
        
    end
        
end


