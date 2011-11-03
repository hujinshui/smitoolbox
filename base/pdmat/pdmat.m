function S = pdmat(a1, a2, a3)
% Constructs a struct representing a (semi)-positive definite matrix
%
%   A pdmat struct can represent a (semi)positive definite matrix
%   in three different forms, each form has a single-char tag to 
%   identify it.
%
%   - 's':  scalar-form:
%           use a scalar to represent a matrix in form of c * eye(d).
%   - 'd':  diagonal-form:
%           use a column-vector of diagonal entries to represent a
%           diagonal matrix.
%   - 'f':  full-form:
%           directly represent the matrix using its full matlab form.
%
%   Note that one can pack multiple matrices of the same size into
%   a single pdmat struct. 
%
%   A pdmat struct is comprised of the following fields:
%
%   - tag:  a string, which equals 'pdmat'
%   - ty:   the char to indicate which form it is in
%   - d:    the dimension of the matrix (the size of matrix is d x d)
%   - n:    the number of matrices packed in the struct
%   - v:    the value in the specified form to represent the matrix
%           - scalar-form: v is a scalar or a 1 x n row vector
%           - diagonal-form: v is a d x 1 vector or a d x n matrix
%           - full-form: v is a d x d matrix or a d x d x n cube
%
%   S = pdmat(v);
%       creates a pdmat struct to represent a single matrix in a proper
%       form.
%
%       - If v is a scalar, then it creates a 1 x 1 matrix in scalar-form.
%       - If v is a d x 1 vector (d > 1), then it creates a d x d matrix 
%         in diagonal-form.
%       - If v is a d x d matrix (d > 1), then it creates a full-form.
%   
%   S = pdmat(ty, d, v);
%       creates a pdmat struct to represent a single matrix or multiple
%       matrices in the specified form.
%
%       Input arguments:
%       - ty:   the type of the representation form: 's', 'd', or 'f'
%       - d:    the dimension of the matrix/matrices
%       - v:    the representation.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% single argument with form inference

if isfloat(a1) && nargin == 1
    
    v = a1;
    
    if ~(isreal(v) && ndims(v) == 2)
        error('pdmat:invalidarg', ...
            'In single-arg input, v should be a real scalar/vector/matrix.');
    end
    
    S.tag = 'pdmat';
    
    if isscalar(v)
        S.ty = 's';
        S.d = 1;
        
    elseif size(v, 2) == 1
        S.ty = 'd';
        S.d = size(v, 1);
        
    elseif size(v, 2) == size(v, 1)
        S.ty = 'f';
        S.d = size(v, 1);
        
    else
        error('pdmat:invalidarg', 'The size of v is invalid.');
    end
        
    S.n = 1;
    S.v = v;
    
        
%% full specification
    
elseif ischar(a1) && nargin == 3
    
    ty = a1;
    d = a2;
    v = a3;
    
    if ~(ischar(ty) && isscalar(ty))
        error('pdmat:invalidarg', 'ty should be a char scalar.');
    end
    
    if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1)
        error('pdmat:invalidarg', 'd should be a positive integer scalar.');
    end
    
    if ~(isreal(v) && ndims(v) <= 3)
        error('pdmat:invalidarg', ...
            'v should be a real scalar/vector/matrix/cube.');
    end
    
    S.tag = 'pdmat';
    
    switch ty
        case 's'
            if ndims(v) == 2 && size(v, 1) == 1
                S.ty = 's';
                S.d = d;
                S.n = size(v, 2);
                S.v = v;
            else
                error('pdmat:invalidarg', ...
                    'For scalar-form, v should be a scalar or a row vector.');
            end                
            
        case 'd'
            if ndims(v) == 2 && size(v, 1) == d
                S.ty = 'd';
                S.d = d;
                S.n = size(v, 2);
                S.v = v;
            else
                error('pdmat:invalidarg', ...
                    'For diagonal-form, the size of v should be d x 1 or d x n.');
            end
            
        case 'f'
            if size(v, 1) == d && size(v, 2) == d
                S.ty = 'f';
                S.d = d;
                S.n = size(v, 3);
                S.v = v;
            else
                error('pdmat:invalidarg', ...
                    'For full-form, the size of v should be d x d or d x d x n.');
            end
            
        otherwise
            error('pdmat:invalidarg', 'ty is invalid.');
    end
    

%% otherwise, error    
    
else
    error('pdmat:invalidarg', 'The input arguments to pdmat are invalid.');
end


