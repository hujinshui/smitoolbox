classdef gaussmrf
    % The class to represent a Gaussian MRF model
    %
    
    % Created by Dahua Lin, on Oct 21, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        nnodes; % the number of nodes
        nedges; % the number of edges        
        graph;  % the underlying graph struct (gr_adjlist)
        
        dims;   % the dimensions of all nodes
        tdim;   % the total dimension
        npots;  % the number of potential vectors of each node
        
        hs;     % the potential sub-vector of all nodes
        Js;     % the information sub-matrix of all nodes
        Rs;     % the information cross-matrix (corresponding to edges)        
    end
    
    methods
        
        function obj = gaussmrf(hs, Js, edges, Rs)
            % Constructs a Gaussian MRF object
            %
            %   obj = gaussmrf(hs, Js, edges, Rs);
            %       constructs a Gauss-MRF model using information-form
            %       parameterization. 
            %
            %       Suppose the MRF has n nodes and m edges. Then
            %       hs and Js should be both cell array with n cells.
            %       In particular, hs{i} and Js{i} are respectively
            %       the potential vector and information matrix
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
            
            % verify hs and Js
            
            if ~(iscell(hs) && isvector(hs))
                error('gaussmrf:invalidarg', 'hs should be a cell vector.');
            end
            n = numel(hs);
            
            if ~(iscell(Js) && isvector(Js) && numel(Js) == n)
                error('gaussmrf:invalidarg', 'Js should be a cell vector with n cells.');
            end
            
            % verify the contents of hs and Js
            
            ds = zeros(n, 1);
            K = size(hs{1}, 2);            
            for i = 1 : n
                h = hs{i};
                J = Js{i};
                
                if ~(isfloat(h) && ndims(h) == 2 && size(h, 2) == K)
                    error('gaussmrf:invalidarg', ...
                        'some potential vectors in hs are invalid.');
                end
                d = size(h, 1);
                
                if ~(isfloat(J) && isequal(size(J), [d d]))
                    error('gaussmrf:invalidarg', ...
                        'some information matrices in Js are invalid.');
                end
                
                ds(i) = d;
            end
           
            % verify graph and Rs
            
            if ~(isnumeric(edges) && ndims(edges) == 2 && size(edges,2) == 2)
                error('gaussmrf:invalidarg', ...
                    'edges should be a matrix of size m x 2.');
            end            
            m = size(edges, 1);
            
            for k = 1 : m
                sd = ds(edges(k, 1));
                td = ds(edges(k, 2));
                
                R = Rs{k};
                if ~(isfloat(R) && isequal(size(R), [sd td]))
                    error('gaussmrf:invalidarg', ...
                        'some matrices in Rs are invalid.');
                end
            end
            
            % build graph and extend Rs
            
            G = gr_adjlist('u', n, edges);
                        
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
            obj.npots = K;
            
            obj.hs = hs;
            obj.Js = Js;
            obj.Rs = Rs;            
        end
        
        
        function h = potential_vector(obj)
            % Get the full potential vector of the model
            %
            %   h = obj.potential_vector;
            %
            
            h = vertcat(obj.hs{:});            
        end
        
        
        function J = information_matrix(obj)
            % Get the full information matrix of the entire model
            %
            %   J = obj.precision_matrix;
            %
            % Note: the resultant matrix is a sparse matrix.
            %
            
            ds = obj.dims;           
            ibase = [0; cumsum(ds(1:end-1))];
            
            n = obj.nnodes;
            srcs = obj.graph.s;
            tars = obj.graph.t;
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
                s = srcs(k) + 1;
                t = tars(k) + 1;                
                [i, j, w] = find(Rs_{k});
                
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

end



