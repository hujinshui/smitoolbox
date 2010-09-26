classdef L2mrf_sfix
    % The class to represent a semi-fixed L2 MRF
    %
    
    % Created by Dahua Lin, on Sep 21, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        mrf0;       % the underlying MRF over all nodes
        nnodes;     % the number of all nodes (including fixed and unfixed)
        
        is_fixed;   % a boolean vector indicating which nodes are fixed
        fixs;       % the indices of fixed nodes
        unfixs;     % the indices of unfixed nodes
                
        Wuu;        % the affinity matrix over unfixed nodes        
        Wuf;        % the relation matrix between fixed and unfixed nodes
        La;         % the augmented Laplacian matrix for inference [n_u x n_u]
        
        Xfix = [];  % the values at fixed nodes (can be [], if unknown)
        Ffix = [];  % the matrix: Xf * Wuf' (can be [], if Xfix is unknown)        
    end
    
    
    methods
        
        function obj = L2mrf_sfix(mrf0, finds, fvals)
            % Constructs a semi-fixed L2 Markov random field
            %
            %   obj = L2mrf_sfix(mrf0, finds);
            %       constructs a semi-fixed L2 Markov random field based 
            %       on mrf0, which is the L2 MRF over all nodes.
            %
            %       The nodes whose indices are given by finds are fixed.
            %
            %   obj = L2mrf_sfix(mrf0, finds, fvals);
            %       with this syntax, one can also specifies the the 
            %       the values at fixed nodes. fvals(:,i) gives the
            %       value at node with index finds(i).
            %
            
            if ~isa(mrf0, 'L2mrf')
                error('L2mrf_sfix:invalidarg', 'mrf0 should be an L2mrf object.');
            end
            
            if ~(isnumeric(finds) && isvector(finds))
                error('L2mrf_sfix:invalidarg', 'finds should be a numeric vector.');
            end
                                                
            if size(finds, 1) > 1  % ensure that finds is a row vector
                finds = finds.';
            end
            
            n = mrf0.nnodes;
            fx = false(1, n);
            fx(finds) = true;
            uinds = find(~fx);
            
            W0 = mrf0.W;
            Wu = W0(uinds, uinds);
            R = W0(uinds, finds);
            Q = laplacemat(Wu, sum(R, 2));
            
            obj.mrf0 = mrf0;
            obj.nnodes = n;
           
            obj.is_fixed = fx;
            obj.fixs = finds;
            obj.unfixs = uinds;
            
            obj.Wuu = Wu;
            obj.Wuf = R;
            obj.La = Q;
            
            if nargin >= 3 && ~isempty(fvals)
                                
                if ~(isfloat(fvals) && ndims(fvals) == 2)
                    error('L2mrf_sfix:invalidarg', 'fvals should be a numeric matrix.');
                end
                
                if size(fvals,2) ~= numel(finds)
                    error('L2mrf_sfix:invalidarg', 'The size of fvals is invalid.');
                end
                
                obj.Xfix = fvals;
                obj.Ffix = fvals * R';
            end
        end
        
        
        function Xu = solve(obj, Xf, Hu, a, fsolver)
            % Solve the values at unfixed nodes
            %            
            %   Xu = obj.solve(Xf, Hu, a);
            %   Xu = obj.solve(Xf, Hu, a, fsolver);
            %       solves the values at unfixed nodes and returns them
            %       as columns in Xu. Xu(:,i) corresponds to obj.unfixs(i).
            %
            %       In the input, Xf gives the values at fixed nodes. 
            %       It can be empty (only when obj.Xfix is given) or a
            %       numeric matrix with nf columns, where nf is the number
            %       of fixed nodes.
            %
            %       Hu gives external forces enforced at each unfixed
            %       nodes (adding a term (- <Hu, Xu>) to the objective
            %       function to be minized). Hu can be empty or omitted,
            %       then, no external forces are enforced.
            %
            %       a gives the regularization coefficient at each
            %       unfixed node (adding a term a_i * ||x_i||^2 / 2).
            %       Here, a can be either a scalar or a vector of 
            %       length nu. a can be empty or omitted.
            %
            %       In addition, the caller can specify the linear equation
            %       solver to use.
            %
            
            % process Xf
            
            if nargin < 2 || isempty(Xf)
                Xf = obj.Xfix;
                if isempty(Xf)
                    error('L2mrf_sfix:solve:invalidarg', ...
                        'Xf must be non-empty when obj.Xfix is not given.');
                end
                Ff = obj.Ffix;
            else
                nf = numel(obj.fixs);
                if ~(isfloat(Xf) && ndims(Xf) == 2)
                    error('L2mrf_sfix:solve:invalidarg', ...
                        'Xf should be a numer matrix.');
                end
                if size(Xf, 2) ~= nf
                    error('L2mrf_sfix:solve:invalidarg', ...
                        'The size of Xf is invalid.');
                end
                Ff = Xf * obj.Wuf';
            end
            
            % process Hu
            
            if nargin < 3 || isempty(Hu)
                H = Ff;
            else
                H = Ff + Hu;
            end
            
            % regularize
            
            Q = obj.La;
            if nargin >= 4 && ~isempty(a)
                nu = numel(obj.unfixs);
                if ~(isfloat(a) && (isscalar(a) || isequal(size(a), [1, nu])))
                    error('L2mrf_sfix:solve:invalidarg', ...
                        'The regularization coefficient a is invalid.');
                end                
                if isscalar(a)
                    if a ~= 0
                        Q = spdiag(nu, a) + Q;
                    end
                else
                    Q = spdiag(a) + Q;
                end
            end
            
            % solve            
            
            if nargin < 5 || isempty(fsolver)
                Xu = H / Q;
            else
                Xu = fsolver(Q, H')';
            end

        end
        
        
    end
    
end