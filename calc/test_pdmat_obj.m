function test_pdmat_obj(obj)
% Test the consistency of a pdmat object
%
%   test_pdmat_obj(obj);
%

% Created by Dahua Lin, on Mar 23, 2010
%

%% verify class

assert(isa(obj, 'pdmat'), 'test_pdmat_obj:invalidarg', ...
    'The input object should be of class pdmat.');

%% Test fullform and is_compatible

d = obj.dim;
m = obj.num;

C = obj.fullform();

% verify size
if m == 1
    assert(isequal(size(C), [d d]));
else
    assert(isequal(size(C), [d d m]));
end

% verify symmetry and positive definiteness
for i = 1 : m
    cC = C(:,:,i);
    assert(isequal(cC, cC'));
    assert(all(eig(cC) >= 0));
end

assert(is_subtype(obj, obj));


%% Test take

for i = 1 : m
    cobj = obj.take(i);
    assert(cobj.dim == d);
    assert(cobj.num == 1);
    
    assert(isequal(cobj.fullform, C(:,:,i)));
end

%% Test inv

invobj = inv(obj);
assert(isa(invobj, class(obj)));
assert(invobj.dim == d);
assert(invobj.num == m);

invC = invobj.fullform;

for i = 1 : m
    devcheck('inv-times-obj', invC(:,:,i) * C(:,:,i), eye(d), 1e-10);
    devcheck('obj-times-inv', C(:,:,i) * invC(:,:,i), eye(d), 1e-10);    
end

%% Test cmult and mtimes and mldivide and cldiv

X = randn(d, m);
Y = cmult(obj, X);

Y0 = zeros(d, m);
for i = 1 : m
    Y0(:,i) = C(:,:,i) * X(:,i);
end

devcheck('cmult', Y, Y0, 1e-12);

n = 20;
X = rand(d, n);

if m == 1
    Y = obj * X;
    Y0 = C * X;        
else
    Y = zeros(d, n, m);
    Y0 = zeros(d, n, m);
    for i = 1 : m
        Y(:,:,i) = obj.take(i) * X;
        Y0(:,:,i) = C(:,:,i) * X;
    end
end

devcheck('mtimes', Y, Y0, 1e-12);

if m == 1
    Z = obj \ X;
    Z0 = C \ X;
else
    Z = zeros(d, n, m);
    Z0 = zeros(d, n, m);
    for i = 1 : m
        Z(:,:,i) = obj.take(i) \ X;
        Z0(:,:,i) = C(:,:,i) \ X;
    end
end

devcheck('mldivide', Z, Z0, 1e-12);

X = randn(d, m);
U = cldiv(obj, X);
U0 = zeros(d, m);
for i = 1 : m
    U0(:,i) = C(:,:,i) \ X(:,i);
end

devcheck('cldiv', U, U0, 1e-12);


%% Test linear operations

c = rand();

M0 = C * c;
o = c * obj;
assert(isa(o, class(obj)));
assert(o.dim == obj.dim && o.num == obj.num);

M = o.fullform;
devcheck('scalar-times (m - 1)', M, M0, 1e-14);

if m == 1
    ca = rand(1, 3);
    oa = ca * obj;
    assert(oa.dim == obj.dim && oa.num == 3);
    Ma = fullform(oa);
        
    Ma0 = zeros(d, d, 3);
    for i = 1 : 3
        Ma0(:,:,i) = C * ca(i);
    end
    
    devcheck('scalar-times (1 - m)', Ma, Ma0, 1e-14);
end

cm = rand(1, m);
om = cm * obj;
assert(isa(om, class(obj)));
assert(om.dim == obj.dim && om.num == obj.num);

Mm0 = zeros(d, d, m);
for i = 1 : m
    Mm0(:,:,i) = cm(i) * C(:,:,i);
end
Mm = om.fullform;

devcheck('scalar-times (m - m)', Mm, Mm0, 1e-14);


A0 = C + M0;
oa = obj + o;
A = fullform(oa);

devcheck('addition', A, A0, 1e-15);


if m == 1
    o = combine(obj, 1);
    assert(isa(o, class(obj)));    
    assert(isequal(o.fullform, obj.fullform));
end

cc = rand(1, m);
o = combine(obj, cc);

assert(isa(o, class(obj)));
assert(o.dim == obj.dim);
assert(o.num == 1);
combC = zeros(d, d);
for i = 1 : m
    combC = combC + cc(i) * C(:,:,i);
end

devcheck('combine', combC, o.fullform, 1e-15);


%% Test quadterm

n = 20;
X = rand(d, n);

M = obj.quadterm(X);
assert(isequal(size(M), [m n]));

M0 = zeros(m, n);
for i = 1 : m
    M0(i, :) = sum(X .* (C(:,:,i) * X), 1);
end

devcheck('quadterm', M, M0, 1e-12);


%% Test chol, sqrtm, trace and logdet

L = chol(obj);
assert(isa(L, 'cholmat'));
assert(L.dim == d);
assert(L.num == m);

S = sqrtm(obj);
assert(isa(S, class(obj)));
assert(S.dim == obj.dim);
assert(S.num == obj.num);

v = logdet(obj);
assert(isequal(size(v), [1 m]));

tr = trace(obj);
assert(isequal(size(tr), [1 m]));

L0 = zeros(d, d, m);
S0 = zeros(d, d, m);
v0 = zeros(1, m);
tr0 = zeros(1, m);

for i = 1 : m
    cC = C(:,:,i);
    cL = chol(cC, 'lower');
    L0(:,:,i) = cL;
    S0(:,:,i) = sqrtm(cC);
    v0(i) = 2 * sum(log(diag(cL)));
    tr0(i) = trace(cC);
end

devcheck('chol', L.fullform, L0, 1e-9);
devcheck('sqrtm', S.fullform, S0, 1e-9);
devcheck('logdet', v, v0, 1e-9);
devcheck('trace', tr, tr0, 1e-14);


%% Auxiliary functions

function devcheck(title, X1, X2, thres)

D = X1 - X2;
dev = max(abs(D(:)));

if dev > thres
    warning('test_pdmat_obj:devcheck', ...
        '%s checking with dev %g', title, dev);
end




