classdef gmrf_blk
    % The class to represent a Gaussian MRF model with block variables
    %
    
    % Created by Dahua Lin, on Oct 21, 2010
    % Modified by Dahua Lin, on Feb 12, 2011
    
    properties(GetAccess='public', SetAccess='private')
        
        nnodes; % the number of block nodes
        nedges; % the number of inter-block edges        
        graph;  % the underlying graph struct (gr_adjlist)
        
        dims;   % the vector of block dimension
        tdim;   % the total dimension
        
        Js;     % the sub information matrix of each block
        Rs;     % the cross information matrix between blocks      
    end
    
    methods
        
        function obj = gmrf_blk(Js, edges, Rs)
            % Constructs a Gaussian MRF object
            %
            %   obj = gmrf_blk(Js, edges, Rs);
            %       constructs a Gauss-MRF model with block variables.
            %
            %       Suppose the MRF has n block nodes and m inter-block 
            %       edges. Then Js should be a cell array with n cells.
            %       In particular, Js{i} is the information matrix
            %       for the i-th node.
            %
            %       edges should be an m x 2 array, each row specifies
            %       an edge between two directly connected nodes. 
            %       Then Rs should be a cell array with m cells, and
            %       Rs{k} corresponds to the k-th edge in edges.
            %       Suppose the k-th edge is (s, t), then Rs{k} should
            %       be a ds x dt matrix. 
            %
            %       Note that for each pair of connected nodes, only
            %       (s, t) or (t, s) need to be given in edges.            
            %       
            
            % verify Js
            
            if ~(iscell(Js) && isvector(Js))
                error('gmrf_blk:invalidarg', 'Js should be a cell vector.');
            end
            n = numel(Js);
                        
            ds = zeros(n, 1);         
            for i = 1 : n
                J = Js{i};                
                d = size(J, 1);                
                if ~(isfloat(J) && isequal(size(J), [d d]))
                    error('gmrf_blk:invalidarg', ...
                        'some information matrices in Js are invalid.');
                end
                
                ds(i) = d;
            end
           
            % verify graph and Rs
            
            if ~(isnumeric(edges) && ndims(edges) == 2 && size(edges,2) == 2)
                error('gmrf_blk:invalidarg', ...
                    'edges should be a matrix of size m x 2.');
            end            
            m = size(edges, 1);
            
            for k = 1 : m
                sd = ds(edges(k, 1));
                td = ds(edges(k, 2));
                
                R = Rs{k};
                if ~(isfloat(R) && isequal(size(R), [sd td]))
                    error('gmrf_blk:invalidarg', ...
                        'some matrices in Rs are invalid.');
                end
            end
            
            % build graph and extend Rs
            
            G = gr_adjlist.from_edges('u', n, edges);
                        
            Rs{2 * m} = [];            
            for k = 1 : m
                Rs{m+k} = Rs{k}';
            end
            
            % set fields
            
            obj.nnodes = n;
            obj.nedges = m;
            obj.graph = G;
            
            obj.dims = ds;
            obj.tdim = sum(ds);
            
            obj.Js = Js;
            obj.Rs = Rs;            
        end
        
               
        function J = information_matrix(obj)
            % Get the full information matrix of the entire model
            %
            %   J = obj.information_matrix;
            %
            % Note: the resultant matrix is a sparse matrix.
            %
            
            ds = obj.dims;           
            ibase = [0; cumsum(ds(1:end-1))];
            
            n = obj.nnodes;
            srcs = obj.graph.source_vs;
            tars = obj.graph.target_vs;
            m = length(srcs);
            
            Js_ = obj.Js;
            Rs_ = obj.Rs;
                        
            i_s = cell(n+m, 1);
            j_s = cell(n+m, 1);
            w_s = cell(n+m, 1);
            
            for k = 1 : n
                [i, j, w] = find(Js_{k});
                
                i_s{k} = i + ibase(k);
                j_s{k} = j + ibase(k);
                w_s{k} = w;                
            end
            
            for k = 1 : m
                s = srcs(k);
                t = tars(k);                
                [i, j, w] = find(Rs_{k});
                if size(i, 2) > 1
                    i = i.';
                    j = j.';
                    w = w.';
                end                
                i_s{n+k} = i + ibase(s);
                j_s{n+k} = j + ibase(t);
                w_s{n+k} = w;
            end
            
            ii = vertcat(i_s{:});
            jj = vertcat(j_s{:});
            ww = vertcat(w_s{:});
            
            td = obj.tdim;            
            J = sparse(ii, jj, ww, td, td);   
        end
        
    end
    
    
    methods(Static)
        
        function gm = from_groupvars(J, grps)
            % Constructs a Gaussian MRF by grouping variables
            %            
            %   gm = gmrf_blk.from_groupvars(J, grps);
            %       constructs a Gaussian MRF model with the original
            %       information matrix J over all variables, and
            %       a grouping of variables into blocks specified by grps.
            %
            %       Suppose we want to divide all variables into 
            %       K blocks, then grps should be a cell array with
            %       K cells, where grps{k} gives the indices of 
            %       variables which are assigned to the k-th group.
            %
            
            % verify variables
            
            n = size(J, 1);
            if ~(isfloat(J) && ndims(J) == 2 && n == size(J, 2))
                error('gmrf_blk:from_groupvars:invalidarg', ...
                    'J should be a square matrix of float type.');
            end
            if ~(iscell(grps) && isvector(grps))
                error('gmrf_blk:from_groupvars:invalidarg', ...
                    'grps should be a cell vector.');
            end            
            K = numel(grps);
            
            % construct skeleton graph
                                
            gmap = zeros(n, 1);
            for k = 1 : K
                gmap(grps{k}) = k;
            end
            
            [i, j] = find(J);            
            gi = gmap(i);
            gj = gmap(j);
            
            sge = find(gi < gj);
            gi = gi(sge);
            gj = gj(sge);
            
            edges = unique([gi gj], 'rows');
            m = size(edges, 1);
            
            % extract corresponding matrices
            
            Js_ = cell(K, 1);
            Rs_ = cell(m, 1);
            
            for i = 1 : K                
                Js_{i} = full(J(grps{i}, grps{i}));
            end
            
            for i = 1 : m
                Rs_{i} = full(J(grps{edges(i, 1)}, grps{edges(i, 2)}));
            end
            
            gm = gaussmrf(Js_, edges, Rs_);
        end
        
    end

end



