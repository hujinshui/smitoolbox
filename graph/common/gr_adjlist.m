classdef gr_adjlist < gr_edgelist
    % The class to represent a graph using adjacency list
    %
    
    % Created by Dahua Lin, on Nov 12, 2010
    %
    
    
    properties(GetAccess='public', SetAccess='protected')
        
        o_ds;   % the out-degree of all nodes [n x 1 uint32]
        o_os;   % the section offset for all nodes [n x 1 int32 zero-based]
        o_es;   % the concatenated outgoing edge array [m'x1 int32 zero-based]
        o_ns;   % the concatenated outgoing neighbor array [m'x1 int32 zero-based]
        
        % for directed graph, m' = m, for undirected graph m' = 2m.
    end
    
    
    methods
        %% Getters
        
        function r = out_degrees(G)
            % Get the outgoing degrees of all nodes
            %
            %   G.out_degrees;
            %
            
            r = G.o_ds;
        end
        
        
        function d = out_degree_of(G, i)
            % Get the outgoing degree of a particular vertex i
            %
            %   d = G.out_degree_of(i);
            %
            
            d = G.o_ds(i);            
        end
        
        
        function r = out_edges_of(G, i)
            % Get the list of outgoing edges from a particular vertex i
            %
            %   G.out_edges_of(i);
            %
            
            b = G.o_os(i);
            d = G.o_ds(i);            
            r = G.o_es(b+(1:d)) + 1;            
        end
        
        
        function [vs, eds] = out_neighbors_of(G, i)
            % Get the list of outgoing neighbors of a particular vertex i
            %
            %   vs = G.out_neighbors_of(i);
            %   [vs, eds] = G.out_neighbors_of(i);
            %
            
            b = G.o_os(i);
            d = G.o_ds(i);
            vs = G.o_ns(b+(1:d)) + 1;
            if nargout >= 2
                eds = G.o_es(b+(1:d)) + 1;
            end            
        end
                
    end
    
    
    methods
        function G = set_weights(G, w)
            % Reset the weights of the edges
            %
            %   G = G.set_weights(w);
            %
            
            if ~(isfloat(w) && isvector(w))
                error('gr_adjlist:invalidarg', ...
                    'w should be a numeric vector.');
            end
                        
            m = G.ne;
            if isscalar(w)
                G.ew = constmat(2*m, 1, w);
            elseif isequal(size(w), [m 1])
                G.ew = [w; w];
            elseif isequal(size(w), [2*m 1])
                G.ew = w;
            else
                error('gr_adjlist:invalidarg', ...
                    'The size of w is invalid.');
            end
        end
        
    end    
    
    
    methods(Static)
        
        function G = from_amat(dty, A)
            % Constructs a graph with adjacency list from affinity matrix
            %
            %   G = gr_adjlist.from_amat(dty, A);
            %
            %       constructs a graph based on the affinity matrix A.
            %       Here dty is a char, which can be 'd' or 'u', 
            %       respectively indicating a directed or undirected
            %       graph.
            %
            
            edlist = gr_edgelist.from_amat(dty, A);
            G = gr_adjlist.from_base(edlist);            
        end
        
        
        function G = from_edges(dty, n, varargin)
            % Constructs a graph with adjacency list from given edges
            %
            %   G = gr_adjlist.from_edges(dty, n, [s, t]);
            %   G = gr_adjlist.from_edges(dty, n, [s, t, w]);
            %   G = gr_adjlist.from_edges(dty, n, s, t);
            %   G = gr_adjlist.from_edges(dty, n, s, t, w);
            %
            %       constructs a graph based on the given edges.
            %       Here, n is the number of vertices, s, t, and w
            %       are respectively the source vertices, target vertices,
            %       and weights of edges.
            %
            
            edlist = gr_edgelist.from_edges(dty, n, varargin{:});
            G = gr_adjlist.from_base(edlist);
        end
        
        
        function G = from_base(edlist)
            % Construct a gr_adjlist object from a gr_edgelist object
            %
            %   G = gr_adjlist.from_base(edlist);
            %
            
            % get info
            
            n = edlist.nv;
            s0 = edlist.es;
            t0 = edlist.et;
            w0 = edlist.ew;
            
            if edlist.is_directed
                s = s0;
                t = t0;
                w = w0;
            else                
                s = [s0; t0];
                t = [t0; s0];
                w = [w0; w0];
            end            
            
            % build neighborhood 
            
            if ~isempty(s)            
                [s_, es] = sort(s);
                ns = t(es);
                
                p = valueseg(s_);
                pds = [diff(p); numel(s_)-p(end)+1];
                
                ds = zeros(n, 1);
                ds(s_(p)+1) = pds;
                
                os = [0; cumsum(ds(1:end-1))];
                
                ds = uint32(ds);
                os = int32(os);
                es = int32(es) - 1;
                ns = int32(ns);
            else
                ds = uint32([]);
                os = int32([]);
                es = int32([]);
                ns = int32([]);
            end
            
            
            % set fields
            
            G = gr_adjlist();
            
            G.dtype = edlist.dtype;            
            G.nv = edlist.nv;
            G.ne = edlist.ne;
            
            G.es = s;
            G.et = t;
            G.ew = w;
            
            G.o_ds = ds;
            G.o_os = os;
            G.o_es = es;
            G.o_ns = ns;                
            
        end
                        
    end
    
    
    
    
    methods        
        
        %% Info dump
        
        function dump(G)
            % dump information about the graph
            %
            %   dump(G);
            %
            
            n = G.nv;
            m = G.ne;
            
            if G.is_directed
                fprintf('Directed Adjacency List: \n');
            else
                fprintf('Undirected Adjacency List: \n');
            end
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
            
            fprintf('  out neighbors:\n');
            for i = 1 : n
                d = G.out_degree_of(i);
                [nvs, nes] = G.out_neighbors_of(i);
                
                fprintf('    [%d] (deg = %d): ', i, d);
                
                for j = 1 : d
                    v = nvs(j);
                    e = nes(j);
                    
                    if e <= m
                        fprintf('%d>%d ', e, v);
                    else
                        fprintf('%d<%d ', e-m, v);
                    end
                end
                fprintf('\n');
            end
            fprintf('\n');
            
        end
                
        
    end
      
    
    
end



