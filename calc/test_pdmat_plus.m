function test_pdmat_plus(cname1, cname2, cfgs)
% Test the correctness of implementation of plus between two pdmat objects
%
%   test_gmat_plus(cname1, cname2, cfgs);
%       Test the correctness of implementation of the plus of cname1 and
%       cname2, which are the names of compatible derived classes of pdmat.
%
%       cfgs are n x 3 matrices, where each row is [d m1 m2].
%       d: the matrix size is d x d
%       m1: the number of matrices contained in the 1st addend
%       m2: the number of matrices contained in the 2nd addend
%       

% Created by Dahua Lin, on Apr 5, 2010
% Modified by Dahua Lin, on Apr 13, 2010
%

%% verify input arguments

assert(ischar(cname1) && ismember('pdmat', superclasses(cname1)), ...
    'test_pdmat_plus:invalidarg', ...
    'The 1st argument should be the name of a derived class of pdmat');

assert(ischar(cname2) && ismember('pdmat', superclasses(cname2)), ...
    'test_pdmat_plus:invalidarg', ...
    'The 2nd argument should be the name of a derived class of pdmat');

assert(isnumeric(cfgs) && ndims(cfgs) == 2 && size(cfgs, 2) == 3, ...
    'test_pdmat:invalidarg', ...
    'The 2nd argument should be a n x 3 numeric matrix.');

%% main

for i = 1 : size(cfgs, 1)
    
    d = cfgs(i, 1);
    m1 = cfgs(i, 2);
    m2 = cfgs(i, 3);
    
    assert(m1 == 1 || m2 == 1 || m1 == m2);
    mc = max(m1, m2);
    
    a = feval([cname1 '.random'], d, m1);
    b = feval([cname2 '.random'], d, m2);
    
    fprintf('Test %s + %s [d = %d, m1 = %d, m2 = %d]\n', ...
        cname1, cname2, d, m1, m2);
    
    assert(is_subtype(a, b) || is_subtype(b, a));    
    c1 = a + b;
    c2 = b + a;
    
    % basic check
    
    assert(isa(c1, class(a)) || isa(c1, class(b)));
    assert(c1.dim == d);
    assert(c1.num == mc);
    
    assert(isa(c2, class(a)) || isa(c2, class(b)));
    assert(c2.dim == d);
    assert(c2.num == mc);    
    
    % quantitative check
    
    MA = fullform(a);
    MB = fullform(b);
    MC1 = fullform(c1);
    MC2 = fullform(c2);
    
    if m1 == m2
        MC0 = MA + MB;
    else
        MC0 = bsxfun(@plus, MA, MB);
    end
    
    dev1 = max(abs(MC1(:) - MC0(:)));
    dev2 = max(abs(MC2(:) - MC0(:)));
    if dev1 > 1e-14 || dev2 > 1e-14;
        warning('test_gmat_plus:largedev', ...
            'Large deviation occurred with dev = %g.', dev);
    end
end

