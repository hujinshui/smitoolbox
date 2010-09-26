classdef subspace_base
%The base class for linear subspace models
%
%   Linear subspace models are widely used in machine learning
%   applications. In this class, the model is characterized by
%   an orthonormal basis matrix, with each column represents
%   a dimension.
%   
%   The model can be utilized as an affine transformer, which
%   transforms the input vectors by projecting them onto the
%   subspace with
%       y = B' * (x - x_anchor)
%
%   Suppose, the input space that x lies in is a d-dimensional
%   space, and the subspace that y is in is a q-dimensional
%   space. Then the basis matrix B should be a d x q orthonormal
%   matrix, such that B' * B == eye(q). x_anchor represents the
%   transform origin in the input space.
%

%   History
%       - Created by Dahua Lin, on May 31, 2008
%

    properties(GetAccess = 'public', SetAccess = 'private')
        basis;              % the basis matrix (a d x q orthonormal matrix)
        
        in_dim;             % the dimensionality of input space (d)
        out_dim;            % the dimensionality of output space (q)
        
        x_anchor;           % the anchor of transform in the input space
        xanchor_at_origin;  % whether the x_anchor is a zero vector.                
    end          
    
    methods
        % construction
        
        function obj = subspace_base(B, xa)
            % constructs the subspace object
            %
            %   obj = subspace_base(B);
            %       constructs the subspace model by specifying the
            %       basis matrix B. The x_anchor is set at origin.
            %
            %       For a model with d-dimensional input space, and
            %       q-dimensional subspace (q <= d), B should be an
            %       orthonormal matrix of size d x q, with each
            %       column giving a base vector of unit length.
            %
            %   obj = subspace_base(B, xa);
            %       constructs the subspace model by specifyin both
            %       the basis matrix B, as well as the anchor in
            %       the input space (xa).
            %
            %       Generally, xa is a d x 1 vector. It can also 
            %       be empty or zero, which indicates that the 
            %       anchor is at origin.
            %
            %   Remarks
            %       - It is the caller's responsibility to ensure
            %         that B is an orthonormal matrix, with B' * B
            %         (approximately) equaling an identity matrix.
            %
            %         For the sake of efficiency, the constructor
            %         does not check this.
            %
            
            assert(ndims(B) == 2, 'subspace_base:invalidarg', ...
                'B should be a matrix.');
            
            [d, q] = size(B);
            
            assert(d >= q, 'subspace_base:invalidarg', ...
                'The dimension of subspace should not exceed that of the input space.');
            
            obj.basis = B;
            
            obj.in_dim = d;
            obj.out_dim = q;                        
            
            if nargin <= 1 || isempty(xa) || isequal(xa, 0)
                obj.x_anchor = zeros(d, 1);
                obj.xanchor_at_origin = true;
            else
                assert(isequal(size(xa), [d 1]), 'subspace_base:invalidarg', ...
                    'xa should be a d x 1 vector.');
                
                obj.x_anchor = xa;
                obj.xanchor_at_origin = all(xa == 0);
            end          
            
        end % constructor
        
        
        function sobj = select(obj, inds)
            % creates a subspace model with selected subset of dimensions
            %
            %   sobj = obj.select(inds);
            %       creates a subspace with subset of dimensions specifed
            %       by the indices in inds.
            %
            
            if obj.xanchor_at_origin
                sobj = subspace_base(obj.basis(:, inds));
            else
                sobj = subspace_base(obj.basis(:, inds), obj.x_anchor);
            end
            
        end
        
    end
    
    
    methods
        % transforms
        
        function Y = transform(obj, X, s)
            % projects the input vectors into the subspace
            %
            %   Y = obj.transform(X);
            %       transforms the input vectors by projecting them into
            %       the subspace as
            %           y = B' * (x - x_anchor)
            %
            %       X should be a d x n matrix with each column being
            %       an input vector, where d == in_dim, and n is the 
            %       number of vectors. And, Y will be a q x n matrix 
            %       with its columns being the output vectors.
            %
            %   Y = obj.transform(X, s);
            %       conducts partial projection, in which the vectors
            %       are projected to a subspace spanned by the subset
            %       of basis indexed by s, as
            %           y = B(:, s)' * (x - x_anchor)
            %       
            %       X should be a d x n matrix as in full projection,
            %       while the output Y will be a q' x n matrix, with
            %       q' equaling the number of basis selected by s.
            %
            %       s can be any form of indexing supported by MATLAB,
            %       including integer array or logical array.
            %
                        
            if ~obj.xanchor_at_origin
                X = bsxfun(@minus, X, obj.x_anchor);
            end
            
            if nargin == 2
                Y = obj.basis' * X;
            else
                Y = obj.basis(:, s)' * X;
            end            
            
        end % transform
        
        
        function X = reconstruct(obj, Y, s)
            % reconstructs the vector in input space
            %
            %   X = obj.reconstruct(Y);
            %       reconstructs the vectors in the input space from
            %       the transformed vectors using backward transform,
            %       which uses the values in the transformed vectors
            %       as coefficient of linear combination of basis, as
            %           x = B * y + x_anchor
            %
            %       Y should be a q x n matrix with each column being
            %       a transformed vector in the subspace, while in the
            %       output, X will be a d x n matrix. Here, d and q
            %       respectively denote in_dim and out_dim, and n is the
            %       number of vectors.
            %
            %   X = obj.reconstruct(Y, s);
            %       reconstructs the vectors in the input space from
            %       part of the transformed vectors. Only a subset of
            %       basis is utilized, as
            %           x = B(:, s) * y + x_anchor
            %
            %       In partial reconstruction, Y should be a q' x n 
            %       matrix, where q' is the number of the dimensions
            %       selected by the indexing array s. While, in the
            %       output, X remains a d x n matrix.
            %
            
            if nargin == 2
                X = obj.basis * Y;
            else
                X = obj.basis(:, s) * Y;
            end
            
            if ~obj.xanchor_at_origin
                X = bsxfun(@plus, X, obj.x_anchor);
            end                
            
        end % reconstruct
        
    end
    
end

