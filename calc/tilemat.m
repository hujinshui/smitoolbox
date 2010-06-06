function R = tilemat(Ms)
% Tile multiple matrices to form a single large one
%
%   R = tilemat(Ms);
%       Creates a big matrix by tiling the matrices given in Ms.
%
%       In the input, Ms can be either of the following form:
%       - an array of size p x q x m x n, where M(:,:,i,j) gives a 
%         sub-matrix.
%       - a cell array of size m x n, where each cell is a matrix of
%         size p x q.
%
%       Then the output R is a matrix of size (p x m) x (q x n), in
%       the following form:
%       
%           M_11, M_12, ..., M_1n
%           M_21, M_22, ..., M_2n
%           ..., ..., ..., ...
%           M_m1, M_m2, ..., M_mn
%
%       Here, M_{ij} is the (i, j)-th sub-matrix.
%

% History
% -------
%   - Created by Dahua Lin, on Apr 4, 2010
%   - Change the error handling, on May 24, 2010
%   - Change to support cell array input additionally, on Jun 6, 2010
%

%% main

if isnumeric(Ms)

    if ndims(Ms) > 4
        error('tilemat:invalidarg', ...
            'When it is a numeric array, Ms should have ndims(Ms) <= 4.');
    end
    
    if ndims(Ms) == 2
        R = Ms;
        return;
    else
        [p, q, m, n] = size(Ms); 
    end

elseif iscell(Ms)
    
    if ndims(Ms) > 2
        error('tilemat:invalidarg', ...
            'When it is a cell array, Ms should have ndims(Ms) == 2.');
    end
    
    if numel(Ms) == 1
        R = Ms{1};   
        return;
    else        
        [p, q] = size(Ms{1});
        [m, n] = size(Ms);
        Ms = [Ms{:}];
    end        
    
end


if m == 1
    R = reshape(Ms, [p, q * n]);

else
    R = reshape(Ms, [p, q * m, n]);
    I = reshape(1:q*m, q, m).';
    R = R(:, I(:), :);
    R = reshape(R, [p * m, q * n]);
end


