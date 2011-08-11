function sm = symat(type, varargin)
% Construct a symmetric matrix
%
%   sm = symat('s', d, v);
%
%       Constructs a symat struct that represents one or multiple
%       diagonal matrices in form of v * eye(d).
%
%       Here, d is the dimension of the sample space, and v can be 
%       either a scalar or a row vector (for multi-matrices case).
%       
%   sm = symat('d', v);
%       
%       Constructs a symat struct that represents one or multiple
%       diagonal matrices in form of diag(v);
%
%       Here, v can be a column vector of size d x 1, or a matrix
%       of size d x n (for constructing n matrices).
%
%   sm = symat('f', v);
%       
%       Constructs a symat struct that represents one or multiple
%       full symmetric matrices.
%
%       Here, v can be a symmetric matrix, or a cube with each page
%       being a symmetric matrix.
%

% Created by Dahua Lin, on Aug 11, 2011
%

switch type 
    case 's'
        sm.ty = 's';
        d = varargin{1};
        v = varargin{2};
        
        if ~(isscalar(d) && d == fix(d) && d > 0)
            error('symat:invalidarg', 'Invalid param d for s-type symat.');
        end
        
        if ~(ndims(v) == 2 && size(v,1) == 1 && isreal(v))
            error('symat:invalidarg', 'Invalid param v for s-type symat.');
        end
            
        sm.d = d;
        sm.n = numel(v);
        sm.v = v;
        
    case 'd'
        sm.ty = 'd';
        v = varargin{1};
        
        if ~(ndims(v) == 2 && isreal(v))
            error('symat:invalidarg', 'Invalid param v for d-type symat.');
        end
       
        [sm.d, sm.n] = size(v);
        sm.v = v;
            
    case 'f'
        sm.ty = 'f';
        v = varargin{1};
        
        d = size(v,1);
        if ~(ndims(v) <= 3 && d == size(v, 2) && isreal(v))
            error('symat:invalidarg', 'Invalid param v for f-type symat.');
        end
        
        sm.d = d;
        sm.n = size(v, 3);
        sm.v = v;
        
    otherwise
        error('symat:invalidarg', 'Invalid type name %s', type);
end


