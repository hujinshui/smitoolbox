classdef gr_edgelist
    % The class to represent a graph with edge list
    %
    
    % Created by Dahua Lin, on Nov 12, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')        
        dty;    % the direction type ('d': directed, 'u': undirected)
        nv;     % the number of vertices
        ne;     % the number of edges
        
        es;     % the source vertices of all edges [ne x 1 int32 zero-based]
        et;     % the target vertices of all edges [ne x 1 int32 zero-based]
        ew;     % the edge weights [empty or m x 1 numeric]        
    end                    
    
    
    methods
        %% Getters
        
        function b = is_directed(G)
            % Get whether the graph is directed
            %
            %   s = G.is_directed;
            %
            
            b = G.dty == 'd';
        end        
        
        function s = source_vs(G)
            % Get the source vertices of all edges
            %
            %   s = G.source_vs;
            %       returns an ne x 1 array of source vertices
            %       [int32 one-based]
            
            s = G.es + 1;            
        end
        
        function t = target_vs(G)
            % Get the target vertices of all edges
            %
            %   t = G.target_vs;
            %       returns an ne x 1 array of target vertices
            %       [int32 one-based]
            %
            
            t = G.et + 1;            
        end
        
        function w = weights(G)
            % Get the weights of all edges
            %
            %   w = G.weights;
            %       returns an ne x 1 array of edge weights (for
            %       weighted graph), or empty (for non-weighted)
            %
            
            w = G.ew;
        end
                
        function E = edges(G)
            % Get all edges
            %
            %   E = G.edges;
            %       returns an ne x 2 array of edges, of which
            %       each row corresponds to an edge.
            %
            
            E = [G.es G.et] + 1;
        end        
        
        function b = is_weighted(G)
            % Tests whether the graph has weighted edges
            %
            %   G.is_weighted;
            %
            
            b = ~isempty(G.ew);
        end

    end
    
        
    methods
    
        %% Constructor
        
        function G = gr_edgelist(dty, varargin)
            % Construct a graph object with edge list
            %
            %   G = gr_edgelist('d', A);
            %   G = gr_edgelist('u', A);
            %       constructs an edge list graph from an adjacency matrix.
            %
            %       The first argument indicates whether the graph is
            %       directed ('d'), or undirected ('u').
            %
            %       A here can be either a logical or numeric array. 
            %       If A is logical, then the function creates an edge 
            %       list with unweighted edges, otherwise it creates an 
            %       edge list with weighted edges.
            %
            %       For undirected graph, only those edges (s, t) with 
            %       s < t are preserved.
            %
            %   G = gr_edgelist(dty, n, [s, t]);
            %   G = gr_edgelist(dty, n, [s, t, w]);
            %   G = gr_edgelist(dty, n, s, t);
            %   G = gr_edgelist(dty, n, s, t, w);
            %       constructs an edge list graph from explicitly given
            %       edges. 
            %
            %       dty here can be either 'd' or 'u' indicating
            %       directed and undirected graph, respectively.
            %
            %       Here, s, t are source and target node indices, and 
            %       w are corresponding edge weights.
            % 
            
            % verify dty
            
            if ~(ischar(dty) && isscalar(dty) && (dty == 'u' || dty == 'd'))
                error('gr_edgelist:invalidarg', 'dty is invalid.');
            end

            % do construction
            
            if nargin == 2 % from affinity matrix
                
                A = varargin{1};                
                n = size(A, 1);
                if ~(isnumeric(A) && ndims(A) == 2 && n == size(A,2))
                    error('gr_edgelist:invalidarg', ...
                        'A should be a square numeric matrix.');
                end
                
                if isnumeric(A)
                    [s, t, w] = find(A);
                else
                    [s, t] = find(A);
                    w = [];
                end
                
                if dty == 'u'
                    se = find(s < t);
                    s = s(se);
                    t = t(se);
                    
                    if ~isempty(w)
                        w = w(se);
                    end
                end
                
                m = numel(s);   
                
            elseif nargin == 3 % from edge array
                
                n = varargin{1};
                E = varargin{2};
                
                if ~(isnumeric(n) && isscalar(n) && n >= 0)
                    error('gr_edgelist:invalidarg', ...
                        'n should be a non-negative scalar.');
                end
                n = double(n);
                
                if ~(ndims(E) == 2 && isnumeric(E))
                    error('gr_edgelist:invalidarg', ...
                        'The 2nd argument should be a numeric matrix.');
                end
                
                m = size(E, 1);
                
                if size(E,2) == 2
                    s = E(:,1);
                    t = E(:,2);
                    w = [];
                    
                elseif size(E,2) == 3
                    s = E(:,1);
                    t = E(:,2);
                    w = E(:,3);
                    
                else
                    error('gr_edgelist:invalidarg', ...
                        'The 2nd argument should have 2 or 3 columns.');
                end
                
            elseif nargin == 4 || nargin == 5 % from vertex lists
                
                n = varargin{1};
                s = varargin{2};
                t = varargin{3};
                
                if nargin < 5
                    w = [];
                else
                    w = varargin{4};
                end
                
                if ~(isnumeric(n) && isscalar(n) && n >= 0)
                    error('gr_edgelist:invalidarg', ...
                        'n should be a non-negative scalar.');
                end
                n = double(n);                
                
                % verify s , t, and w
                if ~(isnumeric(s) && isreal(s) && isvector(s))
                    error('gr_edgelist:invalidarg', 's should be a real vector.');
                end
                
                if ~(isnumeric(t) && isreal(t) && isvector(t))
                    error('gr_edgelist:invalidarg', 't should be a real vector.');
                end
                
                m = numel(s);
                
                if issparse(s); s = full(s); end
                if issparse(t); t = full(t); end
                
                if size(s,2) > 1; s = s.'; end
                if size(t,2) > 1; t = t.'; end
                
                if ~isempty(w)
                    if ~(isnumeric(w) && isreal(w) && isvector(w) && numel(w) == m)
                        error('gr_edgelist:invalidarg', ...
                            'w should be a real vector of length m.');
                    end
                    
                    if issparse(w); w = full(w); end
                    if size(w,2) > 1; w = w.'; end
                end
                
            else
                error('gr_edgelist:invalidarg', ...
                    'The number of input arguments are invalid.');
            end            
            
            G.dty = dty;
            G.nv = n;
            G.ne = m;
            G.es = int32(s) - 1;
            G.et = int32(t) - 1;
            G.ew = w;
                        
        end
        
        
        %% Info dump
        
        function dump(G)
            % dump information about the graph
            %
            %   dump(G);
            %
            
            n = G.nv;
            m = G.ne;
            
            fprintf('Edge List: \n');
            fprintf('------------------------\n');
            fprintf('    # nodes = %d\n', n);
            fprintf('    # edges = %d\n', m);
            fprintf('\n');
                        
            fprintf('  edges: \n');
            s = G.source_vs;
            t = G.target_vs;
            if ~G.is_weighted                
                for i = 1 : m
                    fprintf('    [%d]: (%d, %d)\n', i, s(i), t(i));
                end
            else
                w = G.weights;
                for i = 1 : m
                    fprintf('    [%d]: (%d, %d) = %g\n', i, s(i), t(i), w(i));
                end
            end
            fprintf('\n');
            
        end
        
        
    end
    
    
end


