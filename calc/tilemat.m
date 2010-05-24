function R = tilemat(Ms)
% Tiling the matrices to form a big matrix
%
%   R = tilemat(Ms);
%       Creates a big matrix by tiling the matrices given by pages of 
%       the array Ms.
%   
%       Let Ms be an array of size p x q x m x n, and M_ij denote
%       the matrix given by Ms(:,:,i,j). 
%       Then the output R is a matrix of size (p x m) x (q x n), in
%       the following form:
%       
%           M_11, M_12, ..., M_1n
%           M_21, M_22, ..., M_2n
%           ..., ..., ..., ...
%           M_m1, M_m2, ..., M_mn
%

% Created by Dahua Lin, on Apr 4, 2010
% Change the error handling, on May 24, 2010
%

%% main

if ndims(Ms) > 4
    error('tilemat:invalidarg', ...
        'Ms should be an array with ndims(Ms) <= 4.');
end

if ndims(Ms) == 2
    R = Ms;
    
else
    [p, q, m, n] = size(Ms);
    
    if m == 1
        R = reshape(Ms, [p, q * n]);
        
    else
        R = reshape(Ms, [p, q * m, n]);
        I = reshape(1:q*m, q, m).';
        R = R(:, I(:), :);
        R = reshape(R, [p * m, q * n]);
    end
end


