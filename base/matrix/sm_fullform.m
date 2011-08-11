function M = sm_fullform(sm, i)
% Get the full matrix form of a matrix in an symat struct
%
%   M = sm_fullform(sm, i);
%       Get the full form of the i-th matrix in the symat struct sm.
%
%       If sm is comprised of a single matrix, then i can be omitted.
%       

% Created by Dahua Lin, on Aug 11, 2011
%


% take the specified part of v

if nargin < 2
    if sm.n > 1 
        error('sm_fullform:invalidarg', 'The index of the matrix is needed.');
    end
    v = sm.v;
else
    if ~(isscalar(i) && isnumeric(i))
        error('sm_fullform:invalidarg', 'The index must be a numeric scalar.');
    end
    
    if sm.n == 1 && i == 1
        v = sm.v;
    else    
        v = [];
    end
end

% make it into full-form

switch sm.ty
    case 's'
        if isempty(v)
            v = sm.v(i);
        end
        M = diag(ones(1, sm.d) * v);
    case 'd'
        if isempty(v)
            v = sm.v(:,i);
        end
        M = diag(v);
    case 'f'
        if isempty(v)
            v = sm.v(:,:,i);
        end
        M = v;
end


    
    