classdef tsuite_pdmat
    % The test suite for pdmat_* functions
    %
    
    % Created by Dahua Lin, on Aug 25, 2011
    % Modified by Dahua Lin, on Sep 3, 2011
    
    %% Properties
    
    properties
        types = {'s', 'd', 'f'};
        dims = [1 2 5];
        nums = [1 3];
    end    
        
    %% Test cases
    
    methods
        
        function test_construction(obj)
            
            tys = obj.types;
            ds = obj.dims;
            ns = obj.nums;
            
            % test construction with full spec
            
            for t = 1 : numel(tys)
                ty = tys{t};
                for d = ds
                    for n = ns
                        C = rand_pdmat(ty, d, n, [1 3]);
                        assert(C.d == d);
                        assert(C.n == n);                        
                    end
                end
            end
            
            % test construction with simplified spec
            
            pdm_s = rand_pdmat('s', 1, 1, [1 2]);
            assert(isequal(pdmat(pdm_s.v), pdm_s));
            
            pdm_d = rand_pdmat('d', 2, 1, [1 2]);
            assert(isequal(pdmat(pdm_d.v), pdm_d));
            
            pdm_f = rand_pdmat('f', 2, 1, [1 2]);
            assert(isequal(pdmat(pdm_f.v), pdm_f));
            
        end
        
        function test_pick_and_fullform(obj)
            run_multi(obj, @tsuite_pdmat.do_test_pick_and_fullform);
        end
        
        function test_diag(obj)
            run_multi(obj, @tsuite_pdmat.do_test_diag);
        end
        
        function test_scale(obj)
            run_multi(obj, @tsuite_pdmat.do_test_scale);
        end
        
        function test_inv(obj)
            run_multi(obj, @tsuite_pdmat.do_test_inv);
        end
        
        function test_lndet(obj)
            run_multi(obj, @tsuite_pdmat.do_test_lndet);
        end
        
        function test_mvmul(obj)
            run_multi(obj, @tsuite_pdmat.do_test_mvmul);
        end
                
        function test_lsolve(obj)
            run_multi(obj, @tsuite_pdmat.do_test_lsolve);
        end
        
        function test_quad(obj)
            run_multi(obj, @tsuite_pdmat.do_test_quad);
        end
        
        function test_pwquad(obj)
            run_multi(obj, @tsuite_pdmat.do_test_pwquad);
        end
        
        function test_choltrans(obj)
            run_multi(obj, @tsuite_pdmat.do_test_choltrans);
        end
        
        function test_plus(obj)
            run_multi_pairs(obj, @tsuite_pdmat.do_test_plus);
        end
        
        function test_dot(obj)
            run_multi_pairs(obj, @tsuite_pdmat.do_test_dot);
        end
    end
    
    
    
    methods(Access='private', Static)
        
        function do_test_pick_and_fullform(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            for i = 1 : n
                Si = pdmat_pick(S, i);
                tsuite_pdmat.check_valid(Si, ty, d, 1);
                
                Fi0 = pdmat_fullform(S, i);
                Fi1 = pdmat_fullform(Si);
                assert(isa(Fi0, 'double') && isequal(size(Fi0), [d d]));
                assert(isa(Fi1, 'double') && isequal(size(Fi1), [d d]));
                assert(isequal(Fi0, Fi1));
            end
            
        end
        
        
        function do_test_diag(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            V0 = zeros(d, n);
            for i = 1 : n
                Ci = pdmat_fullform(S, i);
                V0(:, i) = diag(Ci);
            end
            V1 = pdmat_diag(S);
            assert(isequal(V0, V1));
            
        end
        
        
        function do_test_scale(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            c = pi;
            R = pdmat_scale(S, c);
            tsuite_pdmat.check_valid(R, ty, d, n);
            
            for i = 1 : n
                Sm = pdmat_fullform(S, i);
                Rm = pdmat_fullform(R, i);
                assert(norm(Sm * c - Rm) < 1e-15);
            end
            
            if n == 1
                K = 5;
                c = rand(1, K);
                R = pdmat_scale(S, c);
                tsuite_pdmat.check_valid(R, ty, d, K);
                
                Sm = pdmat_fullform(S);
                for i = 1 : K
                    Rm = pdmat_fullform(R, i);
                    assert(norm(Sm * c(i) - Rm) < 1e-15);
                end
                
            else
                c = rand(1, n);
                R = pdmat_scale(S, c);
                tsuite_pdmat.check_valid(R, ty, d, n);
                
                for i = 1 : n
                    Sm = pdmat_fullform(S, i);
                    Rm = pdmat_fullform(R, i);
                    assert(norm(Sm * c(i) - Rm) < 1e-15);
                end
            end
            
        end
        
        
        function do_test_inv(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1 3]);
            R = pdmat_inv(S);
            tsuite_pdmat.check_valid(R, ty, d, n);
            
            for i = 1 : n
                Sm = pdmat_fullform(S, i);
                Rm = pdmat_fullform(R, i);
                
                assert(norm(Sm * Rm - eye(d)) < 1e-12);
            end
            
        end
        
        
        function do_test_lndet(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            r = pdmat_lndet(S);
            assert(isa(r, 'double') && isreal(r));
            assert(isequal(size(r), [1, n]));
            
            for i = 1: n
                Sm = pdmat_fullform(S, i);
                cr0 = lndet(Sm);
                assert(abs(r(i) - cr0) < 1e-12);
            end
        end
        
        
        function do_test_mvmul(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            
            X = rand(d, n);
            Y = pdmat_mvmul(S, X);
            assert(isa(Y, 'double') && isequal(size(Y), [d, n]));
            
            for i = 1 : n
                Sm = pdmat_fullform(S, i);
                y0 = Sm * X(:, i);
                assert(norm(y0 - Y(:,i)) < 1e-13);
            end
            
            if n > 1
                return;
            end
            
            nx = 5;
            X = rand(d, nx);
            Y = pdmat_mvmul(S, X);
            assert(isa(Y, 'double') && isequal(size(Y), [d, nx]));
            
            Y0 = pdmat_fullform(S) * X;
            assert(norm(Y - Y0) < 1e-12);
        end
        
        
        function do_test_lsolve(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            
            X = rand(d, n);
            Y = pdmat_lsolve(S, X);
            assert(isa(Y, 'double') && isequal(size(Y), [d, n]));
            
            for i = 1 : n
                Sm = pdmat_fullform(S, i);
                y0 = Sm \ X(:, i);
                assert(norm(y0 - Y(:,i)) < 1e-12);
            end
            
            if n > 1
                return;
            end
            
            nx = 5;
            X = rand(d, nx);
            Y = pdmat_lsolve(S, X);
            assert(isa(Y, 'double') && isequal(size(Y), [d, nx]));
            
            Y0 = pdmat_fullform(S) \ X;
            assert(norm(Y - Y0) < 1e-12);
        end
        
        
        function do_test_quad(ty, d, n)
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            
            x = rand(d, 1);
            y = rand(d, 1);
            
            rx = pdmat_quad(S, x);
            rxy = pdmat_quad(S, x, y);
            
            assert(isa(rx, 'double') && isequal(size(rx), [n, 1]));
            assert(isa(rxy, 'double') && isequal(size(rxy), [n, 1]));
            
            for i = 1 : n
                Sm = pdmat_fullform(S, i);
                
                crx0 = x' * Sm * x;
                crxy0 = x' * Sm * y;
                
                assert(abs(rx(i) - crx0) < 1e-13);
                assert(abs(rxy(i) - crxy0) < 1e-13);
            end
            
            nx = 5;
            X = rand(d, nx);
            Y = rand(d, nx);
            
            Rxx = pdmat_quad(S, X);
            Rxy = pdmat_quad(S, X, Y);
            
            assert(isa(Rxx, 'double') && isequal(size(Rxx), [n, nx]));
            assert(isa(Rxy, 'double') && isequal(size(Rxy), [n, nx]));
            
            for i = 1 : n
                Sm = pdmat_fullform(S, i);
                
                crx0 = dot(X, Sm * X, 1);
                crxy0 = dot(X, Sm * Y, 1);
                
                assert(norm(Rxx(i, :) - crx0) < 1e-12);
                assert(norm(Rxy(i, :) - crxy0) < 1e-12);
            end
        end
        
        
        function do_test_pwquad(ty, d, n)
            
            if n > 1
                return;
            end
            
            S = rand_pdmat(ty, d, n, [1, 3]);
            
            nx = 5; X = rand(d, nx);
            ny = 6; Y = rand(d, ny);
            
            Rxx = pdmat_pwquad(S, X);
            Rxy = pdmat_pwquad(S, X, Y);
            
            assert(isa(Rxx, 'double') && isequal(size(Rxx), [nx, nx]));
            assert(isa(Rxy, 'double') && isequal(size(Rxy), [nx, ny]));
            
            Sm = pdmat_fullform(S);
            Rxx0 = X' * Sm * X;
            Rxy0 = X' * Sm * Y;
            
            assert(norm(Rxx - Rxx0) < 1e-12);
            assert(norm(Rxy - Rxy0) < 1e-12);
        end
        
        
        function do_test_choltrans(ty, d, n)
            
            if n == 1
                S = rand_pdmat(ty, d, n, [1, 3]);
                
                nx = 5; X = rand(d, nx);
                
                R = pdmat_choltrans(S, X);
                assert(isa(R, 'double') && isequal(size(R), [d, nx]));
                
                Sm = pdmat_fullform(S);
                R0 = chol(Sm, 'lower') * X;
                
                assert(norm(R - R0) < 1e-12);
            end
        end
        

        function do_test_plus(d, ty1, n1, ty2, n2)
            
            S1 = rand_pdmat(ty1, d, n1, [1, 3]);
            S2 = rand_pdmat(ty2, d, n2, [1, 3]);
            
            c1 = 2;
            c2 = 3;
            
            R = pdmat_plus(S1, S2);
            C = pdmat_plus(S1, S2, c1, c2);
            
            nr = max(n1, n2);
            tyr = tsuite_pdmat.max_form_ty(ty1, ty2);
            
            tsuite_pdmat.check_valid(R, tyr, d, nr);
            tsuite_pdmat.check_valid(C, tyr, d, nr);
            
            for i = 1 : nr
                if n1 == 1
                    Sm1 = pdmat_fullform(S1);
                else
                    Sm1 = pdmat_fullform(S1, i);
                end
                
                if n2 == 1
                    Sm2 = pdmat_fullform(S2);
                else
                    Sm2 = pdmat_fullform(S2, i);
                end
                
                Rm = pdmat_fullform(R, i);
                Rm0 = Sm1 + Sm2;
                assert(norm(Rm - Rm0) < 1e-13);
                
                Cm = pdmat_fullform(C, i);
                Cm0 = c1 * Sm1 + c2 * Sm2;
                assert(norm(Cm - Cm0) < 1e-13);
            end
        end
               
        
        function do_test_dot(d, ty1, n1, ty2, n2)
            
            S1 = rand_pdmat(ty1, d, n1, [1, 3]);
            S2 = rand_pdmat(ty2, d, n2, [1, 3]);
            
            V = pdmat_dot(S1, S2);
            n = max(n1, n2);
            
            assert(isa(V, 'double'));
            assert(isequal(size(V), [1, n]));
            
            V0 = zeros(1, n);
            for i = 1 : n
                if n1 == 1
                    Sm1 = pdmat_fullform(S1);
                else
                    Sm1 = pdmat_fullform(S1, i);
                end
                
                if n2 == 1
                    Sm2 = pdmat_fullform(S2);
                else
                    Sm2 = pdmat_fullform(S2, i);
                end
                
                V0(i) = trace(Sm1' * Sm2);
            end
            
            assert(norm(V - V0) < 1e-13);
            
        end
        
    end
    
    
    %% Test running functions
        
    methods
        
        function run_all(obj)
            % Run all tests
            
            test_construction(obj);
            test_diag(obj);            
        end        
                
        function run_multi(obj, tfunc)
            % Run a specific test case under different conditions
                        
            tys = obj.types;
            ds = obj.dims;
            ns = obj.nums;
            
            for t = 1 : numel(tys)
                ty = tys{t};
                for d = ds
                    for n = ns
                        tfunc(ty, d, n);
                    end
                end
            end
        end        
        
        
        function run_multi_pairs(obj, tfunc)
            % Run binary func/ops under different conditions
            
            tys = obj.types;
            ds = obj.dims;
            ns = obj.nums;
                        
            for t1 = 1 : numel(tys)
                for t2 = 1 : numel(tys)
                    
                    ty1 = tys{t1};
                    ty2 = tys{t2};
                    
                    for d = ds
                        for n = ns
                            
                            if n == 1
                                tfunc(d, ty1, 1, ty2, 1);
                            else
                                tfunc(d, ty1, n, ty2, n);
                                tfunc(d, ty1, 1, ty2, n);
                                tfunc(d, ty1, n, ty2, 1);
                            end
                            
                        end
                    end
                end
            end
        end
        
    end
    
    
    %% Auxiliary functions
    
    methods(Static, Access='private')
                
        function check_valid(S, ty, d, n)
            
            assert(is_pdmat(S));
            assert(ischar(S.tag) && strcmp(S.tag, 'pdmat'));
            assert(ischar(S.ty) && strcmp(S.ty, ty));
            assert(isa(S.d, 'double') && isequal(S.d, d));
            assert(isa(S.n, 'double') && isequal(S.n, n));
            
            switch ty
                case 's'
                    assert(isequal(size(S.v), [1, S.n]));
                case 'd'
                    assert(isequal(size(S.v), [S.d, S.n]));
                case 'f'
                    if S.n == 1
                        assert(isequal(size(S.v), [S.d, S.d]));
                    else
                        assert(isequal(size(S.v), [S.d, S.d, S.n]));
                    end
            end
        end
        
        function tyr = max_form_ty(ty1, ty2)
            
            s = [ty1, ty2];
            
            switch s
                case 'ss'
                    tyr = 's';
                case {'sd', 'ds', 'dd'}
                    tyr = 'd';
                case {'sf', 'df', 'fs', 'fd', 'ff'}
                    tyr = 'f';
                otherwise
                    error('test_pdmat:rterror', 'Invalid types.');
            end            
        end
        
    end
    
end

