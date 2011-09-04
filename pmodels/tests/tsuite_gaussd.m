classdef tsuite_gaussd
    % A test suite for gaussd class
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Jun 19, 2010
    %       - Modified by Dahua Lin, on Aug 17, 2011
    %       - Modified by Dahua Lin, on Aug 25, 2011
    %       - Modified by Dahua Lin, on Sep 3, 2011
    %

    %% Properties

    properties
        types = {'s', 'd', 'f'};
        dims = [1 2 5];
        nums = [1, 3];        
    end
    
    %% Test cases
    
    methods
    
        function test_construction(obj)
            run_multi(obj, @tsuite_gaussd.do_test_construction);
        end
                
        function test_statistics(obj)
            run_multi(obj, @tsuite_gaussd.do_test_statistics);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_gaussd.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_gaussd.do_test_sampling);
        end
        
    end
    
    
    %% Test case implementation
    
    methods(Static, Access='private')
        
        function do_test_construction(ty, d, m, is_zeromean, is_shared)
            
            if is_zeromean
                assert(m == 1);
            end
            
            if is_shared
                m2 = 1;
            else
                m2 = m;
            end
            
            if is_zeromean
                mu = zeros(d, m);
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(ty, d, m2, [1, 3]);
            
            g_mp0 = gaussd.from_mp(mu, C);
            g_mp1 = gaussd.from_mp(mu, C, 'ip');
                        
            assert(g_mp0.zmean == is_zeromean);
            assert(g_mp0.shared_cov == (m2 == 1));
            assert(g_mp1.zmean == is_zeromean);
            assert(g_mp1.shared_cov == (m2 == 1));
            
            h = g_mp1.h;
            J = g_mp1.J;
            c0 = g_mp1.c0;
            
            g_cp0 = gaussd.from_ip(h, J, c0);
            g_cp1 = gaussd.from_ip(h, J, c0, 'mp');
            
            g_cp0_a = gaussd.from_ip(h, J);
            g_cp1_a = gaussd.from_ip(h, J, [], 'mp');
                                    
            assert(g_cp0.zmean == is_zeromean);
            assert(g_cp0.shared_cov == (m2 == 1));
            assert(g_cp1.zmean == is_zeromean);
            assert(g_cp1.shared_cov == (m2 == 1));
            
            % basic verification
            
            tsuite_gaussd.basic_verify(g_mp0, d, m, true, false);
            tsuite_gaussd.basic_verify(g_mp1, d, m, true, true);
            tsuite_gaussd.basic_verify(g_cp0, d, m, false, true);
            tsuite_gaussd.basic_verify(g_cp1, d, m, true, true);
            tsuite_gaussd.basic_verify(g_cp0_a, d, m, false, true);
            tsuite_gaussd.basic_verify(g_cp1_a, d, m, true, true);
            
            tsuite_gaussd.verify_mp(g_mp0, mu, C);
            tsuite_gaussd.verify_mp(g_mp1, mu, C);
            tsuite_gaussd.verify_cp(g_cp0, h, J, c0);
            tsuite_gaussd.verify_cp(g_cp1, h, J, c0);
            tsuite_gaussd.verify_cp(g_cp0_a, h, J);
            tsuite_gaussd.verify_cp(g_cp1_a, h, J);
            
            devcheck('c0 consistency', g_cp0_a.c0, c0, 1e-12);
            devcheck('c0 consistency', g_cp1_a.c0, c0, 1e-12);
            
            if is_zeromean
                assert(isequal(g_cp1.mu, 0));
            else            
                devcheck('mu consistency', g_cp1.mu, mu, 1e-12);
            end
            devcheck('C consistency', g_cp1.C.v, C.v, 1e-12);                        
        end
        
        
        function do_test_statistics(ty, d, m, is_zeromean, is_shared)
                                    
            if is_zeromean; assert(m == 1); end            
            if is_shared; m2 = 1; else m2 = m; end            
            if is_zeromean; 
                mu = zeros(d, m); 
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(ty, d, m2, [1, 3]);
            
            g_mp = gaussd.from_mp(mu, C, 'ip');            
            g_cp = gaussd.from_ip(g_mp.h, g_mp.J);

            vars = pdmat_diag(C);
            if is_shared && m > 1
                vars = repmat(vars, [1 m]);
            end
            
            assert(isequal(mean(g_mp), mu));
            assert(isequal(var(g_mp), vars));
            if m == 1 || is_shared
                assert(isequal(cov(g_mp), pdmat_fullform(C)));
            end
            for i = 1 : m
                if is_shared
                    assert(isequal(cov(g_mp, i), pdmat_fullform(C)));
                else
                    assert(isequal(cov(g_mp, i), pdmat_fullform(C, i)));
                end
            end
            
            ents0 = zeros(1, m);
            for i = 1 : m
                if is_shared
                    sig_i = C;
                else
                    sig_i = pdmat_pick(C, i);
                end
                ents0(i) = (1/2) * ( d * log( (2*pi*exp(1)) ) + pdmat_lndet(sig_i) );
            end
            ents1 = entropy(g_mp);
            ents2 = entropy(g_cp);
            
            ents1a = zeros(1, m);
            ents2a = zeros(1, m);
            for i = 1 : m
                ents1a(i) = entropy(g_mp, i);
                ents2a(i) = entropy(g_cp, i);
            end
            
            assert(isequal(size(ents1), [1, m]));
            assert(isequal(size(ents2), [1, m]));
            
            devcheck('entropy calc (C)', ents0, ents1, 1e-12);
            devcheck('entropy calc (J)', ents0, ents2, 1e-12);
            devcheck('entropy calc (C-i)', ents0, ents1a, 1e-12);
            devcheck('entropy calc (J-i)', ents0, ents2a, 1e-12);                        
        end        
        
                
        function do_test_evaluation(ty, d, m, is_zeromean, is_shared)
                        
            if is_zeromean; assert(m == 1); end
            if is_shared; m2 = 1; else m2 = m; end
            if is_zeromean;
                mu = zeros(d, m);
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(ty, d, m2, [1, 3]);            
            g = gaussd.from_mp(mu, C, 'ip');

            if isequal(mu, 0)
                mu = zeros(d, m);
            end
            
            if m == 1
                C2 = pdmat_fullform(g.C);
            else
                C2 = zeros(d, d, m);
                if g.shared_cov
                    sC2 = pdmat_fullform(g.C);
                    C2 = repmat(sC2, [1, 1, m]);
                else
                    for k = 1 : m
                        C2(:,:,k) = pdmat_fullform(g.C, k);
                    end
                end
            end
            
            ld = zeros(1, m);
            Q2 = zeros(d, d, m);
            for k = 1 : m
                ld(k) = lndet(C2(:,:,k));
                Q2(:,:,k) = inv(C2(:,:,k));
            end
            
            n0 = 100;
            X0 = randn(d, n0);
            
            Dsq0 = tsuite_gaussd.comp_sq_mahdist(mu, Q2, X0);
            
            Dsq1 = sqmahdist(g, X0);
            Dsq2 = zeros(m, n0);
            for i = 1 : m
                Dsq2(i, :) = sqmahdist(g, X0, i);
            end
            
            assert(isequal(size(Dsq1), [m n0]));
            devcheck('sqmahdist_1', Dsq0, Dsq1, 1e-8 * max(abs(Dsq0(:))));
            devcheck('sqmahdist_2', Dsq0, Dsq2, 1e-8 * max(abs(Dsq0(:))));
            
            LP0 = bsxfun(@plus, -0.5 * Dsq0, 0.5 * (- ld.' - d * log(2*pi)));
            LP1 = logpdf(g, X0);
            LP2 = zeros(m, n0);
            for i = 1 : m
                LP2(i, :) = logpdf(g, X0, i);
            end
            
            assert(isequal(size(LP1), [m n0]));
            devcheck('logpdf_1', LP0, LP1, 1e-8 * max(abs(LP0(:))));
            devcheck('logpdf_2', LP0, LP2, 1e-8 * max(abs(LP0(:))));
            
            if d == 1
                pint = zeros(1, m);
                for i = 1 : m
                    pint(i) = tsuite_gaussd.prob_integrate_1d(g, i, mu(:,i), Q2(:,:,i));
                end
                devcheck('prob_integrate_1d', pint, ones(1, m), 1e-7);
            elseif d == 2
                pint = zeros(1, m);
                for i = 1 : m
                    pint(i) = tsuite_gaussd.prob_integrate_2d(g, i, mu(:,i), Q2(:,:,i));
                end
                devcheck('prob_integrate_2d', pint, ones(1, m), 1e-5);
            end
        end
        
        
        function do_test_sampling(ty, d, m, is_zeromean, is_shared)
            
            if is_zeromean; assert(m == 1); end
            if is_shared; m2 = 1; else m2 = m; end
            if is_zeromean;
                mu = zeros(d, m);
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(ty, d, m2, [1, 3]);            
            g = gaussd.from_mp(mu, C, 'ip');
            
            C0 = zeros(d, d, m);
            for i = 1 : m
                if is_shared
                    C0(:,:,i) = pdmat_fullform(C);
                else
                    C0(:,:,i) = pdmat_fullform(C, i);
                end
            end
                        
            n = 100000;
            ns = n * ones(1, m);
                                  
            if m == 1
                X = g.sample(ns);
            else
                X = g.sample(ns, 1:m);
            end
            
            mu_s = zeros(d, m);
            C_s = zeros(d, d, m);
            
            for i = 1 : m
                Xi = X(:, (i-1)*n+1 : i*n);
                [C_s(:,:,i), mu_s(:,i)] = veccov(Xi);
            end                        
                        
            devcheck('sample_mean', mu, mu_s, 2e-2);
            devcheck('sample_cov', C0, C_s, 5e-2);
        end
    end    
    
    
    %% Auxiliary functions
    
    methods(Access='private')
        
        function run_multi(obj, tfunc)
            % Run a test case under multiple conditions
            
            tys = obj.types;
            ds = obj.dims;
            ms = obj.nums;
            
            for it = 1 : numel(tys)
                ty = tys{it};
                for d = ds
                    for m = ms
                        tfunc(ty, d, m, false, false);
                        tfunc(ty, d, m, false, true);
                        if m == 1
                            tfunc(ty, d, m, true, false);
                            tfunc(ty, d, m, true, true);
                        end
                    end
                end
            end            
            
        end
    end
    
    
    methods(Static, Access='private')
        
        function basic_verify(g, d, m, has_mp, has_ip)
            
            assert(g.dim == d && g.num == m && ...
                g.has_mp == has_mp && g.has_ip == has_ip);
        end
        
        
        function verify_mp(g, mu, C)
            
            if size(mu,2) == 1 && all(mu == 0)
                assert(isequal(g.mu, 0) && g.zmean)
            else
                assert(isequal(g.mu, mu) && ~g.zmean);
            end
            assert(isequal(g.C, C));
        end
        
        
        function verify_cp(g, h, J, c0)
            
            assert(isequal(g.h, h) && isequal(g.J, J));
            
            if nargin >= 4
                assert(isequal(g.c0, c0));
            end
        end
        
        function D = comp_sq_mahdist(mu, icov, X)
            % compute squared Mahalanobis distance to Gaussian centers
            
            m = size(mu, 2);
            n = size(X, 2);
            
            D = zeros(m, n);
            for i = 1 : m
                A = icov(:,:,i);
                v = mu(:,i);
                D(i, :) = pwsqmahdist(v, X, A);
            end
        end
        
        function v = prob_integrate_1d(g, i, mu, icov)
            
            assert(g.dim == 1);
            s = sqrt(1 / icov);
            
            n = 2000;
            x = linspace(mu - 10 * s, mu + 10 * s, n);
            p = g.pdf(x, i);
            
            v = trapz(x, p);
        end
        
        function v = prob_integrate_2d(g, i, mu, icov)
            
            assert(g.dim == 2);
            C = inv(icov);
            
            [V, D] = eig(C);
            v1 = V(:,1);
            v2 = V(:,2);
            s1 = sqrt(D(1,1));
            s2 = sqrt(D(2,2));
            
            [X1, X2] = meshgrid(-8:0.1:8, -8:0.1:8);
            
            X = bsxfun(@times, X1(:)' * s1, v1) + bsxfun(@times, X2(:)' * s2, v2);
            X = bsxfun(@plus, X, mu);
            
            p = g.pdf(X, i);
            
            v = sum(p) * (s1 * s2) * 0.01;
            
        end
    end
            
end
        

