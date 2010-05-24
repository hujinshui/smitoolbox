function test_gstransform()
% Test the correctness of implementation of gstransform
%

% Created by Dahua Lin, on Apr 6, 2010
%


%% Main

% configuration

gtypes = {'gmat_iso', 'gmat_diag', 'gmat_c2d', 'gmat_comp'};
ds = [1 2 3 5];
qs = [1 2 3 5];
ms = [0 1 3 5];
ns = [1 5 10];
ops = 'NTRC';

cfgs = cartprod(cellfun(@numel, {gtypes, ds, qs, ms, ns, ops}));

cwr = 0;
cwu = 0;
cc = 0;

for i = 1 : size(cfgs, 2)
    c = cfgs(:, i);
    
    gtype = gtypes{c(1)};
    d = ds(c(2));
    q = qs(c(3));
    m = ms(c(4));
    n = ns(c(5));
    op = ops(c(6));
    
    if strcmp(gtype, 'gmat_c2d') && d ~= 2
        continue;
    end
    
    if (op == 'C' && q > 1) || (op == 'R' && d > 1)
        continue;
    end    
    
    cc = cc + 1;
    
    [wr0, wu0] = test_case(gtype, d, q, m, n, op, false);
    [wr1, wu1] = test_case(gtype, d, q, m, n, op, true);
    
    cwr = cwr + wr0 + wr1;
    cwu = cwu + wu0 + wu1;
end

if cwr == 0 && cwu == 0
    fprintf('All %d cases are completed without large deviation warnings.\n', cc);
else
    fprintf('In all %d cases, it raises %d large deviations on R ad %d large deviations on U.\n', cc, cwr, cwu);
end



%% Core test function

function [wr, wu] = test_case(gtype, d, q, m, n, op, tied_A)

fprintf('Test case with g of %s, d = %d, q = %d, m = %d, n = %d, op = %c, tied_A = %d\n', ...
    gtype, d, q, m, n, op, tied_A);

% prepare data

g = feval([gtype '.random'], d, 1);
C = fullform(g);

if m == 0
    W = [];
    m = 1;
else
    W = rand(m, n);
end

if ~tied_A    
    As0 = randn(d, q, n, m);    
    switch op
        case 'N'
            As = As0;
        case 'T'
            As = permute(As0, [2 1 3 4]);
        case 'R'
            As = reshape(As0, [q n m]);
        case 'C'
            As = reshape(As0, [d n m]);
    end
else
    As0 = randn(d, q, n);
    switch op
        case 'N'
            As = As0;
        case 'T'
            As = permute(As0, [2 1 3]);
        case 'R'
            As = reshape(As0, [q n]);
        case 'C'
            As = reshape(As0, [d n]);
    end
    As0 = repmat(As0, [1 1 1 m]);
end

Y = randn(d, n);

% run 

gt0 = gstransform(g, As, W, op);
[gt, U] = gstransform(g, As, Y, W, op);

assert(isa(gt, 'gmat') && gt.num == m && gt.dim == q);
assert(isequal(size(U), [q m]));

R = fullform(gt);

assert(isequal(class(gt0), class(gt)) && isequal(fullform(gt0), R));

% compute ground truth

R0 = zeros(q, q, m);
U0 = zeros(q, m);

for k = 1 : m
    r = 0;
    u = 0;
    
    if isempty(W)
        w = ones(1, n);
    else
        w = W(k, :);
    end
    
    for i = 1 : n
        Aik = As0(:,:,i,k);
        r = r + w(i) * Aik' * C * Aik;
        u = u + w(i) * Aik' * C * Y(:,i);
    end
    
    R0(:,:,k) = r;
    U0(:,k) = u;
end

% quantitative comparison

dev_r = max(abs(R0(:) - R(:)));
dev_u = max(abs(U0(:) - U(:)));

wr = 0;
wu = 0;

if dev_r > 1e-12
    warning('test_gstransform:largedev', ...
        'Large deviation of R with dev_r = %g.', dev_r);
    wr = 1;
end

if dev_u > 1e-12
    warning('test_gstransform:largedev', ...
        'Large deviation of U with dev_u = %g.', dev_u);
    wu = 1;
end









