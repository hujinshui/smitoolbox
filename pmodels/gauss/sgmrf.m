classdef sgmrf
    % The class to represent a simplified Gaussian MRF as follows
    %
    %   The quadratic part of the energy function is given by
    %
    %       E_q(x) 
    %       = (1/2) * 
    %         (
    %           sum_i a_i x_i^2 + 
    %           sum_e w_e (x_{e_i} - x_{e_j})^2 
    %         )
    %       = (1/2) * x^T J x
    %
    %       Here, J is an n x n positive definite matrix, with
    %       J_{ii} = a_i + sum_j w_{ij}, and J_{ij} = -w_{ij}
    %       for i ~= j.
    %       
    
    % Created by Dahua Lin, on Nov 2, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        nnodes;     % the number of nodes (each node is a scalar variable)
        nedges;     % the number of edges
        graph;      % the underlying graph (gr_adjlist)
        
        a;          % the regularization coefficients (n x 1 or scalar)
        J;          % the information matrix (n x n sparse)
    end
        
    methods
        
        function obj = sgmrf(g, a)
            % constructs a simplified Gaussian MRF model
            %
            %   obj = sgmrf(W, a);
            %   obj = sgmrf(g, a);
            %       constructs a simplified Gaussian MRF model.
            %       
            %       In the input, we can specify the underlying 
            %       graph by the affinity matrix W, or a weighted
            %       graph struct g (gr_edgelist or gr_adjlist).
            %
            %       The input argument a is a vector of length n,
            %       where a(i) is the additional regularization
            %       coefficient for the i-th variable.
            %
            %   Remarks
            %   -------
            %       - a should be a non-negative vector, and in
            %         each connected component of the underlying
            %         graph, at least one variable should be
            %         a positive a-value.
            %
            
            % verify input arguments
            
            g = gr_adjlist(g, 'u');
            if isempty(g.w)
                error('sgmrf:invalidarg', 'g should be a weighted graph.');
            end                        
            
            % compute J
            
            J_ = laplacemat(g, a); 
            if size(a, 2) > 1; a = a.'; end
            
            % set fields
            
            obj.nnodes = g.n;
            obj.nedges = g.m;
            obj.graph = g;
            
            obj.a = a;
            obj.J = J_;            
        end
        
        
        function mu = infer(obj, h)
            % infer the mean vector from a potential vector
            %
            %   mu = obj.infer(h);
            %       infers the mean vector mu from the potential 
            %       vector h, by solving J^{-1} h.
            %
            %       h can be a column vector of length n or
            %       multiple vectors organized into an n x K
            %       matrix. In output, x has the same size as y.
            
            mu = obj.J \ h;
        end
                
        function x = smooth(obj, y)
            % Performs smoothing.
            %
            %   x = obj.smooth(y);
            %       it solves the following problem:
            %       
            %           minimize 
            %           (1/2) sum_i a_i (x_i - y_i)^2 + 
            %           (1/2) sum_e w_e (x_{e_i} - x_{e_j})^2.
            %
            %       y can be a column vector of length n or 
            %       multiple vectors organized into an n x K matrix.
            %       In output, x has the same size as y.
            %
            
            n = obj.nnodes;
            if ~(isfloat(y) && ndims(y) == 2 && size(y,1) == n)
                error('sgmrf:smooth:invalidarg', ...
                    'The size of y is invalid.');
            end
            
            a_ = obj.a;                        
            if size(y,2) == 1 || isscalar(a_)
                ay = a_ .* y;
            else
                ay = bsxfun(@times, a_, y);
            end
            
            x = obj.J \ ay;
        end
        
    end
    
end

