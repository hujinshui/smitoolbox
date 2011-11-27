function W = wishart_sample(S, m, n)
% Samples from a Wishart distribution
%
%   W = wishart_sample(d, m);
%       draws a sample from a d-dimensional m-degree standard Wishart
%       distribution (whose scale matrix is an identity matrix).
%
%       The result is returned in form of a d x d matrix.
%
%   W = wishart_sample(d, m, n);
%       draws n samples from a d-dimensional m-degree standard Wishart
%       distribution.
%
%       The results are returned as an d x d x n array, with each page
%       being a sample.
%
%   W = wishart_sample(S, m);
%   W = wishart_sample(S, m, n);
%       draws n samples from an m-degree Wishart distribution with
%       scale matrix S, which should be a pdmat struct (S.n == 1).
%       If n is omitted, it draws one sample
%
%   Remarks
%   -------
%       - The implementation is based on Odell and Feiveson's method,
%         and Bartlett decomposition, and its complexity is O(d^2) for
%         each sample.
%

%   History
%   -------
%       - Created by Dahua Lin, on Spe 2, 2011
%

%% verify input arguments

if isnumeric(S) && isscalar(S)
    d = S;
    if ~(d == fix(d) && d >= 1)
        error('wishart_sample:invalidarg', ...
            'd should be a positive integer scalar.');
    end
    S = [];
    
elseif is_pdmat(S)
    if S.n ~= 1
        error('wishart_sample:invalidarg', ...
            'S should contain only one matrix.');
    end
    d = S.d;
    
else
    error('wishart_sample:invalidarg', ...
        'The first argument to wishart_sample is invalid.');
end

if ~(isfloat(m) && isreal(m) && isscalar(m) && m >= d - 1)
    error('wishart_sample:invalidarg', ...
        'The deg of freedom m should be a real scalar with m >= d - 1.');
end

if nargin < 3
    n = 1;
else
    if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
        error('wishart_sample:invalidarg', ...
            'n should be a non-negative integer scalar.');
    end
end

%% main

if d == 1
    W = randg(0.5 * m, 1, n) * 2;
    
    if ~isempty(S)
        W = S.v * W;
    end
    
    if n ~= 1
        W = reshape(W, [1 1 n]);
    end
    
elseif d == 2
    N = randn(1, n);
    
    V1 = randg(0.5 * m, 1, n) * 2;
    V2 = randg(0.5 * (m-1), 1, n) * 2;
    
    if n == 1
        W = zeros(2, 2);
        W(1) = V1;
        W(2) = N * sqrt(V1);
        W(3) = W(2);
        W(4) = V2 + N.^2;
    else
        W = zeros(4, n);
        W(1, :) = V1;
        W(2, :) = N .* sqrt(V1);
        W(3, :) = W(2, :);
        W(4, :) = V2 + N.^2;
        W = reshape(W, [2, 2, n]);
    end
    
    if ~isempty(S)
        switch S.ty
            case 's'
                W = S.v * W;
            case 'd'
                sqv = sqrt(S.v);
                M = sqv * sqv';
                if n == 1
                    W = M .* W;
                else
                    W = bsxfun(@times, M, W);
                end
            case 'f'
                L = chol2x2(S.v);
                
                if n == 1
                    W = L * W * L';
                    W = (W + W') * 0.5;
                else
                    W = reshape(W, [2, 2 * n]);
                    U = L * W;
                    U = [ ...
                        reshape(U(:, 1:2:end), 2*n, 1), ...
                        reshape(U(:, 2:2:end), 2*n, 1) ];
                    U = U * L';
                    
                    W = zeros(2, 2*n);
                    W(:, 1:2:end) = reshape(U(:, 1), [2, n]);
                    W(:, 2:2:end) = reshape(U(:, 2), [2, n]);
                    W = reshape(W, [2, 2, n]);
                    
                    w12 = (W(1,2,:) + W(2,1,:)) * 0.5;
                    W(1,2,:) = w12;
                    W(2,1,:) = w12;
                end
        end
    end
    
    
elseif d == 3
    
    N = randn(3, n);
    N12 = N(1, :);
    N13 = N(2, :);
    N23 = N(3, :);
    
    V1 = randg(0.5 * m, 1, n) * 2;
    V2 = randg(0.5 * (m-1), 1, n) * 2;
    V3 = randg(0.5 * (m-2), 1, n) * 2;
    
    if n == 1
        W = zeros(3, 3);
        
        W(1) = V1;                      % (1, 1)
        W(5) = V2 + N12.^2;             % (2, 2)
        W(9) = V3 + N13.^2 + N23.^2;    % (3, 3)
        
        W(2) = N12 .* sqrt(V1);         % (2, 1)
        W(4) = W(2);                    % (1, 2)
        
        W(3) = N13 .* sqrt(V1);         % (3, 1)
        W(7) = W(3);                    % (1, 3)
        
        W(6) = N23 .* sqrt(V2) + N12 .* N13;    % (3, 2)
        W(8) = W(6);                            % (2, 3)
    else
        W = zeros(9, n);
        
        W(1,:) = V1;                        % (1, 1)
        W(5,:) = V2 + N12.^2;               % (2, 2)
        W(9,:) = V3 + N13.^2 + N23.^2;      % (3, 3)
        
        W(2,:) = N12 .* sqrt(V1);           % (2, 1)
        W(4,:) = W(2,:);                    % (1, 2)
        
        W(3,:) = N13 .* sqrt(V1);           % (3, 1)
        W(7,:) = W(3,:);                    % (1, 3)
        
        W(6,:) = N23 .* sqrt(V2) + N12 .* N13;      % (3, 2)
        W(8,:) = W(6,:);                            % (2, 3)
        
        W = reshape(W, [3 3 n]);
    end
    
    if ~isempty(S)
        switch S.ty
            case 's'
                W = S.v * W;
            case 'd'
                sqv = sqrt(S.v);
                M = sqv * sqv';
                if n == 1
                    W = M .* W;
                else
                    W = bsxfun(@times, M, W);
                end
            case 'f'
                L = chol(S.v, 'lower');
                
                if n == 1
                    W = L * W * L';
                    W = (W + W') * 0.5;
                else
                    W = reshape(W, [3, 3 * n]);
                    U = L * W;
                    U = [ ...
                        reshape(U(:, 1:3:end), 3*n, 1), ...
                        reshape(U(:, 2:3:end), 3*n, 1), ...
                        reshape(U(:, 3:3:end), 3*n, 1)];
                    U = U * L';
                    
                    W = zeros(3, 3*n);
                    W(:, 1:3:end) = reshape(U(:, 1), [3, n]);
                    W(:, 2:3:end) = reshape(U(:, 2), [3, n]);
                    W(:, 3:3:end) = reshape(U(:, 3), [3, n]);
                    W = reshape(W, [3, 3, n]);
                    
                    w12 = (W(1,2,:) + W(2,1,:)) * 0.5;
                    w13 = (W(1,3,:) + W(3,1,:)) * 0.5;
                    w23 = (W(2,3,:) + W(3,2,:)) * 0.5;
                    
                    W(1,2,:) = w12; W(2,1,:) = w12;
                    W(1,3,:) = w13; W(3,1,:) = w13;
                    W(2,3,:) = w23; W(3,2,:) = w23;
                    
                end
        end
    end
    
else % higher d
    
    % generate N ~ normal
    N = randn(d * (d-1) / 2, n);
    
    % generate V ~ chi-square -> sqrt
    A = 0.5 * ((m + 1) - (1:d)');
    if n ~= 1
        A = A(:, ones(1, n));
    end
    V = sqrt(randg(A) * 2);
    
    % make W
    ltri = tril(true(d, d), -1);
    di = 1:(d+1):(d*d);
    T = zeros(d, d);
    
    if n == 1
        T(ltri) = N;
        T(di) = V;
        
        if ~isempty(S)
            T = pdmat_choltrans(S, T);
        end
        W = T * T';
        
    else
        W = zeros(d, d, n);
        ltri = find(ltri);
        
        if isempty(S)
            for i = 1 : n
                T(ltri) = N(:,i);
                T(di) = V(:, i);
                W(:,:,i) = T * T';
            end
        else
            switch S.ty
                case 's'
                    sqv = sqrt(S.v);
                    for i = 1 : n
                        T(ltri) = N(:,i);
                        T(di) = V(:, i);
                        T = sqv * T;
                        W(:,:,i) = T * T';
                    end
                case 'd'
                    sqv = sqrt(S.v);
                    for i = 1 : n
                        T(ltri) = N(:,i);
                        T(di) = V(:, i);
                        T = bsxfun(@times, sqv, T);
                        W(:,:,i) = T * T';
                    end
                case 'f'
                    L = chol(S.v, 'lower');
                    for i = 1 : n
                        T(ltri) = N(:,i);
                        T(di) = V(:, i);
                        T = L * T;
                        W(:,:,i) = T * T';
                    end
            end
            
        end                
    end
    
end



