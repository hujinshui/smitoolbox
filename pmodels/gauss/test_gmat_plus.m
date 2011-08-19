function test_gmat_plus()
% Test the correctness of the implementation of gmat_plus
%
%   test_gmat_plus;
%

% Created by Dahua Lin, on Aug 18, 2011
%

%% main

cforms = 'sdf';
ds = [1, 2, 5];
ns = [1, 3];

for cf1 = cforms
for cf2 = cforms
    for d = ds
        for n1 = ns
        for n2 = ns
            
            fprintf('Testing on [cf = %c, d = %d, n1 = %d] with [cf = %c, d = %d, n2 = %d]\n', ...
                cf1, d, n1, cf2, d, n2);
            
            C1 = randpdm(cf1, d, n1);
            C2 = randpdm(cf2, d, n2);
            
            F1 = gmat_fullform(cf1, d, C1);
            F2 = gmat_fullform(cf2, d, C2);
            
            [Cr, cfr] = gmat_plus(C1, cf1, C2, cf2);
            
            Fr0 = bsxfun(@plus, F1, F2);
            Fr = gmat_fullform(cfr, d, Cr);
            
            assert(isequal(Fr, Fr0));
            
        end
        end
    end
end
end



%% Auxiliary functions

function C = randpdm(cf, d, m)

switch cf
    case 's'
        C = 0.5 + rand(1, m);
    case 'd'
        C = 0.5 + rand(d, m);
    case 'f'
        C = zeros(d, d, m);
        for i = 1 : m
            R = orth(randn(d));
            dv = 0.5 + rand(d, 1);
            C(:,:,i) = R' * diag(dv) * R;
        end
end


