function L = afpclus(S)
%AFPCLUS Clustering by Affinity Propagation
%
%   L = AFPCLUS(S);
%

L = afpclus_naive(S, 20);


%% core function

function L = afpclus_naive(S, niters)

n = size(S, 1);
A = zeros(n, n);
R = zeros(n, n);

lam = 0.5;

for t = 1 : niters
    
    fprintf('t = %d\n',t);
    
    % compute r
    
    R_pre = R;
    AS = A + S;
    [Y, I] = max(AS, [], 2);
    
    for i = 1 : n
        AS(i, I(i)) = -inf;
    end
    Y2 = max(AS, [], 2);
    
    R = S - Y(:, ones(1, n));
    
    for i = 1 : n
        R(i, I(i)) = S(i, I(i)) - Y2(i);
    end
    
    R = (1 - lam) * R + lam * R_pre;
    
    display(R);
    
    % compute a
    
    A_pre = A;
    Rp = max(R, 0);
    
    for k = 1 : n
        Rp(k, k) = R(k, k);
    end
    
    A = repmat(sum(Rp, 1), [n, 1]) - Rp;
    dA = diag(A);
    A = min(A, 0);
    for k = 1 : n
        A(k, k) = dA(k);
    end
    
    A = (1 - lam) * A + lam * A_pre;    
    
    display(A);
end

E = R + A;
I = find(diag(E) > 0);
K = numel(I);
[~, c] = max(S(:,I), [], 2);
c(I) = 1:K;
L = I(c);
