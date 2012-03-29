function C = mmsum(As, Bs, W, op)
% Compute the (weighted) sum of matrix products
%
%   C = mmsum(As, Bs, W);
%   C = mmsum(As, Bs, W, op);
%       compute the (weighted sum) of products between the matrices in 
%       As and Bs, according to the following formula.
%
%       When op = 'NN':
%           C = \sum_{i=1}^n  w(i) * A_i * B_i,
%       When op = 'NT':
%           C = \sum_{i=1}^n  w(i) * A_i * B_i',
%       When op = 'TN':
%           C = \sum_{i=1}^n  w(i) * A_i' * B_i,
%       When op = 'TT':
%           C = \sum_{i=1}^n  w(i) * A_i' * B_i.
%
%       Here, A_i is given by As(:,:,i) and B_i is given by Bs(:,:,i).
%       The weights are given by each row of W. W can be a matrix of
%       size m x n, then in the output, C has size(C, 3) == 3, and
%       C(:,:,k) is computed based on W(k,:).
%
%       If W is empty, then there is only one group of weights, in 
%       which w(i) always equals 1.      
%
%       The user can input Bs as an empty array to indicate that As and
%       Bs are the same. In such cases, more efficient implementation
%       may be employed.
%
%       If w is omitted or empty, then it assumes w(i) = 1 for all i.
%       If op is omitted, then by default, it is set to 'NN'.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 5, 2010
%       - Modified by Dahua Lin, on Apr 6, 2010
%           - support multiple groups of weights.
%       - Modified by Dahua Lin, on Apr 15, 2010
%           - change error handling.
% 

%% verify input arguments

if nargin < 4
    op = 'NN';
else
    if ~(ischar(op) && numel(op) == 2 && ...
        (op(1) == 'N' || op(1) == 'T') && (op(2) == 'N' || op(2) == 'T'))
        error('mmsum:invalidarg', ...
            'The argument op should be either of ''NN'', ''NT'', ''TN'', ''NN''.');
    end
end

if ~(isfloat(As) && ndims(As) <= 3) 
    error('mmsum:invalidarg', ...
        'As should be a numeric array with ndims(As) <= 3.');
end

if ~isempty(Bs)
    if ~(isfloat(Bs) && ndims(Bs) <= 3) 
        error('mmsum:invalidarg', ...
            'Bs should be either empty or a numeric array with ndims(Bs) <= 3.');
    end

    if size(As, 3) ~= size(Bs, 3)
        error('mmsum:invalidarg', ...
            'As and Bs should have size(As, 3) == size(Bs, 3).');
    end
end

% check the consistency of inner dimension

if op(1) == 'N'
    di_a = size(As, 2);       
else
    di_a = size(As, 1);
end

if op(2) == 'N'
    if isempty(Bs)
        di_b = size(As, 1);
    else
        di_b = size(Bs, 1);
    end
else
    if isempty(Bs)
        di_b = size(As, 2);
    else
        di_b = size(Bs, 2);
    end
end

if di_a ~= di_b
    error('mmsum:invalidarg', ...
        'The inner dimension of As and Bs are inconsistent.');
end

n = size(As, 3);

if ~isempty(W)
    if ~(isfloat(W) && ndims(W) == 2 && size(W,2) == n)
        error('mmsum:invalidarg', 'W should be a numeric matrix with n columns.');
    end
    m = size(W, 1);
else
    m = 1;
end

%% main

if isempty(Bs) && op(1) ~= op(2)  % special case with faster implementation    
    if op(1) == 'N'
        HA = horzc(As);
        C = make_carr(HA, HA, W, n, di_a, 'NT');
    else
        VA = vertc(As);        
        C = make_carr(VA, VA, W, n, di_a, 'TN');
    end
    
    % ensure symmetry
    if m == 1
        C = 0.5 * (C + C');
    else
        for k = 1 : m
            cC = C(:,:,k);
            C(:,:,k) = 0.5 * (cC + cC');
        end
    end
else
    if isempty(Bs)
        Bs = As;
    end
    
    if op(1) == 'N'
        if op(2) == 'N' % op = 'NN'            
            C = make_carr(horzc(As), vertc(Bs), W, n, di_a, 'NN');
        else            % op = 'NT'
            C = make_carr(horzc(As), horzc(Bs), W, n, di_a, 'NT');
        end
    else
        if op(2) == 'N' % op = 'TN'
            C = make_carr(vertc(As), vertc(Bs), W, n, di_a, 'TN');
        else            % op = 'TT'
            C = make_carr(vertc(As), horzc(Bs), W, n, di_a, 'TT');
        end
    end
end
        

%% Auxiliary functions

function Ae = horzc(As)

[p, q, n] = size(As);
if n == 1
    Ae = As;
else
    Ae = reshape(As, p, q * n);
end


function Ae = vertc(As)

[p, q, n] = size(As);
if n == 1
    Ae = As;
else
    I = reshape(1:q*n, q, n).';
    Ae = reshape(As(:, I(:)), p * n, q);
end


function C = make_carr(A, B, W, n, di, op)

if isempty(W)
    C = do_mm(A, B, op);    
else
    if di > 1
        I = reshape(ones(di, 1) * (1:n), 1, n * di);
        W = W(:, I);
    end
    m = size(W, 1);
    
    if op(1) == 'T'
        W = W.';
    end
            
    if m == 1
        C = do_mm(bsxfun(@times, A, W), B, op);
    else    
        if op(1) == 'N'
            C1 = do_mm(bsxfun(@times, A, W(1,:)), B, op);
            C = zeros([size(C1) m], class(C1));
            C(:,:,1) = C1;
            for k = 2 : m
                C(:,:,k) = do_mm(bsxfun(@times, A, W(k,:)), B, op);
            end
        else
            C1 = do_mm(bsxfun(@times, A, W(:, 1)), B, op);
            C = zeros([size(C1) m], class(C1));
            C(:,:,1) = C1;
            for k = 2 : m
                C(:,:,k) = do_mm(bsxfun(@times, A, W(:, k)), B, op);
            end
        end
    end
end


function C = do_mm(A, B, op)

if op(1) == 'N'
    if op(2) == 'N'
        C = A * B;
    else
        C = A * B';
    end
else
    if op(2) == 'N'
        C = A' * B;
    else
        C = (B * A)';
    end
end

      