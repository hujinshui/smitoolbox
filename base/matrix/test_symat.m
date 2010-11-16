function test_symat(classname, ds, ns)
% Test the implementation of symmetric classes
%
%   test_symat(classname);
%   test_symat(classname, ds, ns);
%
%   Test the implementation of a symmetric matrix class.
%

%% verify input

if nargin < 2
    ds = [1 2 5];
end

if nargin < 2
    ns = [1 3 10];
end


%% Main skeleton

test_construction(classname, ds, ns);

test_linear_calc(classname, ds, ns);

test_multiply(classname, ds, ns);

test_left_div(classname, ds, ns);

test_combine(classname, ds, ns);

test_join(classname, ds, ns);

test_inv(classname, ds, ns);

test_character(classname, ds, ns);

test_quad(classname, ds, ns);

test_qtrans(classname, ds);

test_choltrans(classname, ds);


%% core testing functions

function test_construction(classname, ds, ns)
% Test matrix generation, fullform, take, and getm

disp('Test matrix construction and retrieval ...');

for d = ds
    for n = ns
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        A = randpdmat(classname, d, n);
        assert(A.d == d && A.n == n);
        
        M0 = fullform(A);
        assert((n == 1 && isequal(size(M0), [d d])) || ...
            (n > 1 && isequal(size(M0), [d d n])));
        
        for i = 1 : n
            Ai = take(A, i);
            assert(isa(Ai, classname) && Ai.d == d && Ai.n == 1);
            assert(isequal(M0(:,:,i), fullform(Ai), A.getm(i)));
        end        
    end
end


function test_linear_calc(classname, ds, ns)
% Test matrix plus, minus, and scalar multiplication

disp('Test matrix addition and subtraction ...');


for d = ds
    for n = ns
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        k = rand();
        km = rand(1, n);
        A = randpdmat(classname, d, n);
        B = randpdmat(classname, d, n);
        
        S = A + B;
        M = A - B;
        Cl = k .* A;
        Cr = A .* k;
        Clm = km .* A;
        Crm = A .* km;
        
        assert(isa(S, classname) && S.d == d && S.n == n);
        assert(isa(M, classname) && M.d == d && M.n == n);
        assert(isa(Cl, classname) && Cl.d == d && Cl.n == n);  %#ok<BDSCI>
        assert(isa(Cr, classname) && Cr.d == d && Cr.n == n);  %#ok<BDSCI>
        assert(isa(Clm, classname) && Clm.d == d && Clm.n == n);  %#ok<BDSCI>
        assert(isa(Crm, classname) && Crm.d == d && Crm.n == n);  %#ok<BDSCI>
        
        Af = fullform(A);
        Bf = fullform(B);
        Sf = fullform(S);
        Mf = fullform(M);
        Clf = fullform(Cl);
        Crf = fullform(Cr);
        Clmf = fullform(Clm);
        Crmf = fullform(Crm);
        
        checkdev('matrix plus', Af + Bf, Sf, 1e-15);
        checkdev('matrix minus', Af - Bf, Mf, 1e-15);
        
        assert(isequal(Clf, Crf));
        assert(isequal(Clmf, Crmf));
        
        checkdev('matrix scalar multiply', Af * k, Clf, 1e-15);
        checkdev('matrix scalar multiply (m)', ...
            bsxfun(@times, Af, reshape(km, [1 1 n])), Clmf, 1e-15);
    end
end


function test_multiply(classname, ds, ns)
% Test matrix multiplication

disp('Test matrix multiplication ...');

m = 8;
for d = ds
    for n = ns
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        A = randpdmat(classname, d, n);
        B = randpdmat(classname, d, n);
        
        C = A * B;
        assert(isa(C, classname) && C.d == d && C.n == n);
        
        if n == 1        
            L = rand(m, d);
            R = rand(d, m);
            
            LA = L * A;
            AR = A * R;
            
            assert(isnumeric(LA) && isequal(size(LA), [m d]));
            assert(isnumeric(AR) && isequal(size(AR), [d m]));
            
            checkdev('left multiply', L * fullform(A), LA, 1e-12);
            checkdev('right multiply', fullform(A) * R, AR, 1e-12);
        end
        
        X = rand(d, n);
        Y = cmv(A, X);
        assert(isnumeric(Y) && isequal(size(Y), [d n]));
        
        Y0 = zeros(d, n);
        for i = 1 : n
            Y0(:,i) = A.getm(i) * X(:,i);
        end
        
        checkdev('cmv', Y0, Y, 1e-12);
    end
end


function test_left_div(classname, ds, ns)
% Test left division

disp('Test matrix left division ...');

m = 8;
for d = ds
    for n = ns
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        A = randpdmat(classname, d, n);
        X = rand(d, n);
        Y = cdv(A, X);
        
        assert(isnumeric(Y) && isequal(size(Y), [d n]));
        Y0 = zeros(d, n);
        for i = 1 : n
            Y0(:, i) = A.getm(i) \ X(:,i);
        end
        checkdev('cdv', Y0, Y, 1e-12);
        
        if n == 1
            X = rand(d, m);
            Y = A \ X;
            assert(isnumeric(Y) && isequal(size(Y), [d m]));
            
            checkdev('left division', fullform(A) \ X, Y, 1e-12);
        end        
    end
end


function test_combine(classname, ds, ns)
% Test matrix combination

disp('Test matrix combination');

for d = ds
    for n = ns
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        A = randpdmat(classname, d, n);
        
        w = rand(1, n);
        
        C = combine(A);
        Cw = combine(A, w);
        
        assert(isa(C, classname) && C.d == d && C.n == 1);
        assert(isa(Cw, classname) && Cw.d == d && Cw.n == 1);
        
        C0 = zeros(d, d);
        Cw0 = zeros(d, d);
        for i = 1 : n
            C0 = C0 + A.getm(i);
            Cw0 = Cw0 + A.getm(i) * w(i);
        end
        
        checkdev('combine', fullform(C), C0, 1e-14);
        checkdev('combine (w)', fullform(Cw), Cw0, 1e-14);                
    end
end


function test_inv(classname, ds, ns)
% Test the matrix inverse

disp('Test matrix inverse ...');

for d = ds
    for n = ns
        fprintf('\tfor d = %d, n= %d\n', d, n);
        
        A = randpdmat(classname, d, n);
        B = inv(A);
        
        assert(isa(B, classname) && B.d == d && B.n == n);
        
        IAB = A * B; %#ok<MINV>
        IBA = B * A; %#ok<MINV>
        
        I = repmat(eye(d), [1 1 n]);
        
        checkdev('inv (I - AB)', fullform(IAB), I, 1e-12);
        checkdev('inv (I - BA)', fullform(IBA), I, 1e-12);        
    end
end




function test_join(classname, ds, ns)
% Test matrix object joining

m = 4;
disp('Test matrix joining ...');

for d = ds
    for n = ns
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        As = cell(1, m);
        for i = 1 : m
            As{i} = randpdmat(classname, d, n);                        
        end
        
        C1 = join(As{1});
        Cm = join(As{:});
        
        C10 = fullform(As{1});
        G = cellfun(@fullform, As, 'UniformOutput', false);
        Cm0 = cat(3, G{:});
        
        assert(isequal(fullform(C1), C10));
        assert(isequal(fullform(Cm), Cm0));
    end
end


function test_character(classname, ds, ns)
% Test lndet and trace
%

disp('Test lndet and trace ...');


for d = ds
    for n = ns
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        A = randpdmat(classname, d, n);
        
        trv = trace(A);
        ldv = lndet(A);
        
        assert(isnumeric(trv) && isequal(size(trv), [1 n]));
        assert(isnumeric(ldv) && isequal(size(ldv), [1 n]));
        
        trv0 = zeros(1, n);
        ldv0 = zeros(1, n);
        
        for i = 1 : n
            X = A.getm(i);
            trv0(i) = trace(X);
            ldv0(i) = log(det(X));
        end
        
        checkdev('trace', trv, trv0, 1e-15);
        checkdev('lndet', ldv, ldv0, 1e-12);        
    end
end



function test_quad(classname, ds, ns)
% Test the computation of quadratic terms

disp('Test quad ...');

m = 8;

for d = ds
    for n = ns        
        fprintf('\tfor d = %d, n = %d\n', d, n);
        
        A = randpdmat(classname, d, n);        
        X = rand(d, m);
        Y = rand(d, m);
        
        Q = quad(A, X, Y);
        assert(isnumeric(Q) && isequal(size(Q), [n, m]));
        
        Q0 = zeros(n, m);
        for i = 1 : n
            cA = A.getm(i);            
            Q0(i, :) = sum(X .* (cA * Y), 1);            
        end        
        
        checkdev('quad', Q, Q0, 1e-13);        
    end
end


function test_qtrans(classname, ds)
% Test the computation of quadratic transform

disp('Test qtrans ...');

p = 8;

n = 1;
for d = ds       
    fprintf('\tfor d = %d\n', d);
    
    A = randpdmat(classname, d, n);
    B = rand(p, d);
    
    R = qtrans(A, B);
    assert(isnumeric(R) && isequal(size(R), [p, p]));
    
    R0 = B * fullform(A) * B';
    
    checkdev('qtrans', R, R0, 1e-13);
end



function test_choltrans(classname, ds)

disp('Test choltrans ...');

for d = ds
    
    fprintf('\tfor d = %d\n', d);
    
    A = randpdmat(classname, d, 1);
    L0 = chol(fullform(A), 'lower');
    L = A.choltrans(eye(d));
    
    checkdev('choltrans', L, L0, 1e-12);
end

        
%% Auxiliary functions

function A = randpdmat(classname, d, n)

A = feval([classname '.randpdm'], d, n);


function checkdev(item, X, Y, thres)

dev = max(abs(X(:) - Y(:)));
if dev > thres
    warning('test_symat:checkdev', ...
        'Large deviation occurred for %s with dev = %g.\n', item, dev);
end






