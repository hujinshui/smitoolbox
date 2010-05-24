function test_pdmat(cname, cfgs)
% Tests the correctness of implementation of a derived class of pdmat
%
%   test_pdmat(cname, cfgs);
%       tests the correctness of the implementation of cname, which
%       is the name of a derived class of pdmat.
%
%       cfgs is a n x 2 matrix, with each row giving a pair of d and m,
%       which indicates the testing condition (dim == d and num == m).
%
%       The test verifies the consistency of the implementation in
%       various ways.
%
%   Remarks
%       - In order to use this function for testing, the class cname
%         should implement a static method called random that
%         support the following syntax to generate random objects:
%
%           R = cname.random(d, m);
%
%         d is the dimension, and m is the number of matrices contained
%         in the output object R.
%

% Created by Dahua Lin, on Mar 23, 2010
%

%% verify input arguments

assert(ischar(cname) && ismember('pdmat', superclasses(cname)), ...
    'test_pdmat:invalidarg', ...
    'The 1st argument should be the name of a derived class of pdmat');

assert(isnumeric(cfgs) && ndims(cfgs) == 2 && size(cfgs, 2) == 2, ...
    'test_pdmat:invalidarg', ...
    'The 2nd argument should be a n x 2 numeric matrix.');

%% main

for i = 1 : size(cfgs, 1)
    
    d = cfgs(i, 1);
    m = cfgs(i, 2);
    
    obj = feval([cname '.random'], d, m);
    
    fprintf('Test %s [d = %d, m = %d]\n', cname, d, m);
    
    assert(obj.dim == d);
    assert(obj.num == m);
    
    test_pdmat_obj(obj);
end

