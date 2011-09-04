classdef gr_tree
    % The class to represent a tree (in graph theoretical sense)
    %
    
    % Created by Dahua Lin, on Nov 13, 2010
    %
   
    
    properties(GetAccess='public', SetAccess='protected')
        
        nv;     % the number of vertices
        ne;     % the number of edges (nv - ntr)
        ntr;    % the number of connected trees
        
        rts;    % the roots [1 x ntr int32 zero-based]
        pas;    % the parents [n x 1 int32 zero-based]
        
        chds;   % the concatenated array of children [ne x 1 int32 zero-based]
        nchs;   % the number of children of each node [n x 1 int32]
        ofs;    % the offsets of each node in chds [n x 1 int32 zero-based]        
    end
    
    
    methods
        
        function rs = roots(T)
            % Get the root vertices
            %
            %   T.roots;
            %
            
            rs = T.rts + 1;            
        end
        
        
        function ps = parents(T)
            % Get the array of all parents
            %
            %   T.parents;
            %
            
            ps = T.pas + 1;
        end
                
        
        function p = parent_of(T, i)
            % Get the parent of a particular vertex i
            %
            %   T.parent_of(i);
            %
            
            p = T.pas(i) + 1;
        end
        
        
        function cs = children_of(T, i)
            % Get the children vertices of a particular vertex i
            %
            %   T.children_of(i);
            %
            
            cs = T.chds(T.ofs(i)+(1:T.nchs(i))) + 1;
        end
        
        
        function b = is_root(T, i)
            % Test whether a vertex is a root vertex
            %
            %   T.is_root(i);
            %
            
            b = T.pas(i) < 0;
        end
        
        
        function b = is_leaf(T, i)
            % Test whether a vertex is a leaf vertex
            %
            %   T.is_leaf(i);
            %
            
            b = T.nchs(i) == 0;
        end
                        
    end
    
    
    methods(Static)
       
        function T = from_parents(pas)
            % Construct a tree from parents of all vertices
            %
            %   T = gr_tree.from_parents(parents);
            %
            
            if ~(isnumeric(pas) && isvector(pas))
                error('gr_tree:invalidarg', 'parents should be a vector.');
            end
            
            n = numel(pas);
            
            if issparse(pas); pas = full(pas); end
            if size(pas, 2) > 1; pas = pas.'; end
            
            rs = find(pas == 0);
            
            if isempty(rs)
                error('gr_tree:rterror', 'No root is found in the tree.');
            end
            nt = numel(rs);
            
            [sp, cs] = sort(pas); %#ok<ASGLU>
            cs = cs(nt+1:end);
            
            ncs = intcount([1, n], pas).';
            os = [0; cumsum(ncs(1:end-1))];
            
            % set fields
            
            T = gr_tree();
            
            T.nv = n;
            T.ne = n - nt;
            T.ntr = nt;
            
            T.rts = int32(rs) - 1;
            T.pas = int32(pas) - 1;
            
            T.chds = int32(cs) - 1;
            T.nchs = int32(ncs);
            T.ofs = int32(os);
            
        end
                                
    end
    
    
    methods
       
        function dump(T)
            % Dump the tree
            %
            %   dump(T);
            %
            
            n = T.nv;
            
            fprintf('Tree \n');
            fprintf('------------------------\n');
            fprintf('    # nodes = %d\n', n);
            fprintf('    # comps = %d\n', T.ntr);
            fprintf('    # edges = %d\n', T.ne);            
            fprintf('\n');
            
            fprintf('  nodes:\n');
            
            for i = 1 : n
                nc = T.nchs(i); 
                pa = T.parent_of(i);
                fprintf('    [%d] (pa = %d)', i, pa);
                
                if nc > 0
                    fprintf(' ==> ');
                    cs = T.children_of(i);
                    for j = 1 : nc                        
                        fprintf('%d ', cs(j));
                    end
                end
                fprintf('\n');
            end
            fprintf('\n');
            
        end
        
        
    end
    
    
    
end








