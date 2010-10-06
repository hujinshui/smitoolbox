function R = repmrf(W, rw, n)
% Replicate MRFs to introduce inter-copy links
%
%   R = repmrf(W, rw, n);
%       creates a new mrf by replicating n copies of the source mrf whose
%       affinity weights are given by W.  Links between corresponding 
%       notes in different copies are specified by W.
%
%       Each markov random field is represented by a sparse affinity
%       matrix. In the input, W represents the source MRF. 
%       rw specifies the inter-copy links, and n is the number of
%       copies to make.
%
%       Suppose the original MRF has m nodes, then W is a sparse matrix
%       of size m x m, where W(i, j) specifies half of the weights between
%       node i and node j. W(i,i) should be zero for each i.
%
%       rw can be in either of the following forms:
%       - a scalar. Then a link with weight rw is introduced to connect
%         corresponding nodes in adjacent copies.
%       - a vector of length m. Then a link weight rw(i) is introduced to
%         connect between i-th node in adjacent copies of the MRF.
%       - a m x m matrix. Then a link with weight rw(i, j) is introduced 
%         to connect between the i-th node in k-th copy to the j-th
%         node in (k+1)-th copy for each k = 1, ..., n-1.
%         Different from W, rw need not to be symmetry.
%       - empty. No inter-copy links are added.
%

% Created by Dahua Lin, on Apr 16, 2010
%

%% verify and process input arguments

assert(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,1) == size(W,2), ...
    'repmrf:invalidarg', ...
    'W should be a real symmetric matrix.');

m = size(W, 1);

if ~isempty(rw)
    assert(isfloat(rw) && isreal(rw) && ndims(rw) == 2 && ...
        (isscalar(rw) || (isvector(rw) && numel(rw) == m) || isequal(size(rw), [m m])), ...
        'repmrf:invalidarg', 'The argument rw is invalid.');
    
    has_inter = nnz(rw) > 0;
    
    if isscalar(rw)
        rw = rw(ones(m, 1));
    elseif isvector(rw) && size(rw, 2) > 1
        rw = rw.';
    end        
else
    has_inter = false;
end

assert(isscalar(n) && isnumeric(n) && n == fix(n) && n >= 0, ...
    'repmrf:invalidarg', 'n should be a non-negative integer scalar.');

N = m * n; % number of nodes in output

%% main

if n == 0
    R = [];
    
elseif n == 1
    R = W;
    
else % n > 1
    
    % make intra-copy links
    
    [I0, J0, V0] = find(W); 
            
    if n == 2        
        I_intra = [I0; I0 + m];
        J_intra = [J0; J0 + m];
        V_intra = [V0; V0];
                
    else % n > 2
        I_intra = bsxfun(@plus, I0, (0:n-1) * m);
        J_intra = bsxfun(@plus, J0, (0:n-1) * m);
        
        I_intra = I_intra(:);
        J_intra = J_intra(:);
        V_intra = repmat(V0, n, 1);
    end
    
    % make inter-copy links
            
    if has_inter
        
        if isscalar(rw) || isvector(rw)

            I_inter = (1 : m * (n-1)).';
            J_inter = I_inter + m;
            if n == 2
                V_inter = rw;
            else
                V_inter = repmat(rw, n-1, 1);
            end
                        
        else
            
            [I1, J1, V1] = find(rw);
            
            if n == 2
                I_inter = I1;
                J_inter = J1 + m;
                V_inter = V1;
                
            else
                I_inter = bsxfun(@plus, I1, (0:n-2)*m);
                J_inter = bsxfun(@plus, J1, (1:n-1) * m);
                I_inter = I_inter(:);
                J_inter = J_inter(:);
                
                V_inter = repmat(V1, n-1, 1);                
            end                            
        end
        
    end
    
    % combine
            
    if has_inter
        I_a = [I_intra; I_inter; J_inter]; 
        J_a = [J_intra; J_inter; I_inter];
        V_a = [V_intra; V_inter; V_inter];
    else
        I_a = I_intra;
        J_a = J_intra;
        V_a = V_intra;
    end
    
    R = sparse(I_a, J_a, V_a, N, N);
end
    
    