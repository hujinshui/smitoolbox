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
                
        function test_calculation(obj)
            run_multi(obj, @tsuite_gaussd.do_test_calculation);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_gaussd.do_test_evaluation);
        end
        
        function test_mle(obj)
            run_multi(obj, @tsuite_gaussd.do_test_mle);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_gaussd.do_test_sampling);
        end
        
    end
    
    
    %% Test case implementation
    
    methods(Static, Access='private')
        
        function do_test_construction(cf, d, m, is_zeromean, is_shared)
            
            % parse settings
            
            if is_zeromean; assert(m == 1); end            
            if is_shared; m2 = 1; else m2 = m; end
            
            % generate mean-param model
            
            if is_zeromean
                mu = 0;
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(cf, d, m2, [1, 3]);
            
            g_m0 = gaussd('m', mu, C);
            tsuite_gaussd.verify_mp(g_m0, d, m, cf, mu, C);
            
            % generate canon-param model
                        
            if is_zeromean
                h = 0;
            else
                h = pdmat_lsolve(C, mu);
            end
            J = pdmat_inv(C);
                        
            g_c0 = gaussd('c', h, J);
            tsuite_gaussd.verify_cp(g_c0, d, m, cf, h, J);
            
            % mean to canon conversion
                        
            g_c1 = gaussd('c', g_m0);
            
            tsuite_gaussd.verify_cp(g_c1, d, m, cf, [], []);            
            devcheck('mp->cp conversion (h)', g_c1.h, h, 1e-12);
            devcheck('mp->cp conversion (J)', g_c1.J.v, J.v, 1e-12);
            
            % canon to mean conversion
            
            g_m1 = gaussd('m', g_c0);
            
            tsuite_gaussd.verify_mp(g_m1, d, m, cf, [], []);
            devcheck('cp->mp conversion (mu)', g_m1.mu, mu, 1e-12);
            devcheck('cp->mp conversion (C)', g_m1.C.v, C.v, 1e-12);            
            
        end
        
        
        function do_test_calculation(ty, d, m, is_zeromean, is_shared)
                                    
            % parse settings
            
            if is_zeromean; assert(m == 1); end            
            if is_shared; m2 = 1; else m2 = m; end            
            
            % generate model 
            
            if is_zeromean; 
                mu = 0; 
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(ty, d, m2, [1, 3]);
            
            g_mp = gaussd('m', mu, C);            
            g_cp = gaussd('c', g_mp);

            % calculate ground-truths
            
            if m == 1
                if is_zeromean
                    gt.ca = 0;                
                else
                    gt.ca = tsuite_gaussd.calc_ca(mu, C);
                end
            else
                gt.ca = zeros(1, m);
                for k = 1 : m
                    if C.n == 1
                        Ck = C;
                    else
                        Ck = pdmat_pick(C, k);
                    end
                    gt.ca(k) = tsuite_gaussd.calc_ca(mu(:,k), Ck);
                end
            end
            
            if m2 == 1
                gt.cb = tsuite_gaussd.calc_cb(C);
                gt.ent = tsuite_gaussd.calc_entropy(C);
            else
                gt.cb = zeros(1, m);
                gt.ent = zeros(1, m);
                for k = 1 : m
                    Ck = pdmat_pick(C, k);
                    gt.cb(k) = tsuite_gaussd.calc_cb(Ck);
                    gt.ent(k) = tsuite_gaussd.calc_entropy(Ck);
                end
            end
            
            % test ca and cb
            
            [ca_m, cb_m] = gaussd_const(g_mp);
            [ca_c, cb_c] = gaussd_const(g_cp);
            
            devcheck('ca_m', ca_m, gt.ca, 1e-12);
            devcheck('cb_m', cb_m, gt.cb, 1e-12);                                    
            devcheck('ca_c', ca_c, gt.ca, 1e-12);
            devcheck('cb_c', cb_c, gt.cb, 1e-12);  
            
            % test entropy
            
            ent_m = gaussd_entropy(g_mp);
            ent_c = gaussd_entropy(g_cp);
            
            devcheck('ent_m', ent_m, gt.ent, 1e-12);
            devcheck('ent_c', ent_c, gt.ent, 1e-12);            
        end        
        
                
        function do_test_evaluation(ty, d, m, is_zeromean, is_shared)
                        
            % parse settings
            
            if is_zeromean; assert(m == 1); end
            if is_shared; m2 = 1; else m2 = m; end
            
            % generate models
            
            if is_zeromean;
                mu = 0;
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(ty, d, m2, [1, 3]);    
            
            g_mp = gaussd('m', mu, C);
            g_cp = gaussd('c', g_mp);

            % calculate ground-truths            
            
            % generate full form of mu and C
            
            if isequal(mu, 0)
                mu = zeros(d, m);
            end
            
            if m == 1
                C2 = pdmat_fullform(g_mp.C);
            else
                C2 = zeros(d, d, m);
                if is_shared
                    C2 = pdmat_fullform(g_mp.C);
                    C2 = repmat(C2, [1, 1, m]);
                else
                    for k = 1 : m
                        C2(:,:,k) = pdmat_fullform(g_mp.C, k);
                    end
                end
            end
            
            ld = zeros(1, m);
            Q2 = zeros(d, d, m);
            for k = 1 : m
                ld(k) = lndet(C2(:,:,k));
                Q2(:,:,k) = inv(C2(:,:,k));
            end
            
            % generate sample points
            
            n0 = 100;
            X = randn(d, n0);
            
            % evaluate ground-truths
            
            gt.D = tsuite_gaussd.comp_sq_mahdist(mu, Q2, X);    
            ldv = pdmat_lndet(C);
            gt.LP = bsxfun(@plus, -0.5 * gt.D, 0.5 * (- ldv.' - d * log(2*pi)));
            
            % compare with ground-truths
            
            D_m = gaussd_sqmahdist(g_mp, X);
            D_c = gaussd_sqmahdist(g_cp, X);
            
            LP_m = gaussd_logpdf(g_mp, X);
            LP_c = gaussd_logpdf(g_cp, X);
            
            devcheck('sqmahdist (m)', D_m, gt.D, 1e-8 * max(abs(gt.D(:))));
            devcheck('sqmahdist (c)', D_c, gt.D, 1e-8 * max(abs(gt.D(:))));
            
            devcheck('logpdf (m)', LP_m, gt.LP, 1e-8 * max(abs(gt.LP(:))));
            devcheck('logpdf (c)', LP_c, gt.LP, 1e-8 * max(abs(gt.LP(:))));          
        end
        
        
        function do_test_mle(ty, d, m, is_zeromean, is_shared)
            
            % parse settings
            
            if is_zeromean
                return;
            end
            
            % generate data
            
            mu0 = randn(d, 1);
            L0 = rand(d, d);            
            n = 1000;
            X = bsxfun(@plus, mu0, L0 * randn(d, n));
            w = rand(m, n);
            wsp = l2mat(m, randi(m, 1, n), 'sparse');
            
            % perform estimation
            
            if m == 1
                Ge0 = gaussd_mle(X, [], ty, is_shared);                
                tsuite_gaussd.verify_mp(Ge0, d, m, ty, [], []);
                                                
                Ge = gaussd_mle(X, w, ty, is_shared);
                tsuite_gaussd.verify_mp(Ge, d, m, ty, [], []);
                
                Ge2 = gaussd_mle(X, wsp, ty, is_shared);
                tsuite_gaussd.verify_mp(Ge2, d, m, ty, [], []);
                                
                Gr0 = tsuite_gaussd.gmle(X, ones(1, n), ty, is_shared);
                Gr = tsuite_gaussd.gmle(X, w, ty, is_shared);  
                Gr2 = tsuite_gaussd.gmle(X, full(wsp), ty, is_shared);
                
                devcheck('mle (mean)', Ge0.mu, Gr0.mu, 1e-12);
                devcheck('mle (cov)', Ge0.C.v, Gr0.C.v, 1e-12);
                devcheck('w-mle (mean)', Ge.mu, Gr.mu, 1e-12);
                devcheck('w-mle (cov)', Ge.C.v, Gr.C.v, 1e-12);
                devcheck('wsp-mle (mean)', Ge2.mu, Gr2.mu, 1e-12);
                devcheck('wsp-mle (cov)', Ge2.C.v, Gr2.C.v, 1e-12);

            else
                Ge = gaussd_mle(X, w, ty, is_shared);
                tsuite_gaussd.verify_mp(Ge, d, m, ty, [], []);
                
                Ge2 = gaussd_mle(X, wsp, ty, is_shared);
                tsuite_gaussd.verify_mp(Ge2, d, m, ty, [], []);
                
                Gr = tsuite_gaussd.gmle(X, w, ty, is_shared); 
                Gr2 = tsuite_gaussd.gmle(X, full(wsp), ty, is_shared);
                
                devcheck('w-mle (mean)', Ge.mu, Gr.mu, 1e-12);
                devcheck('w-mle (cov)', Ge.C.v, Gr.C.v, 1e-12);
                devcheck('wsp-mle (mean)', Ge2.mu, Gr2.mu, 1e-12);
                devcheck('wsp-mle (cov)', Ge2.C.v, Gr2.C.v, 1e-12);
            end            
            
        end
        
        
        
        function do_test_sampling(ty, d, m, is_zeromean, is_shared)
                        
            % parse settings
                        
            if m > 1
                return;
            end            
            
            if is_zeromean; assert(m == 1); end
            if is_shared; m2 = 1; else m2 = m; end
            
            % generate models
            
            if is_zeromean;
                mu = 0;
            else
                mu = randn(d, m);
            end
            
            C = rand_pdmat(ty, d, m2, [1, 3]);            
            g_mp = gaussd('m', mu, C);
            g_cp = gaussd('c', g_mp);
            
            if is_zeromean
                mu0 = zeros(d, 1);
            else
                mu0 = mu;
            end                            
            C0 = pdmat_fullform(C);
            
            % perform sampling
                        
            n = 100000;
            
            X_m = gaussd_sample(g_mp, n);
            X_c = gaussd_sample(g_cp, n);                                  
                                           
            devcheck('sample_mean (m)', vecmean(X_m), mu0, 2e-2);
            devcheck('sample_cov (m)', veccov(X_m), C0, 5e-2);                        
            devcheck('sample_mean (c)', vecmean(X_c), mu0, 2e-2);
            devcheck('sample_cov (c)', veccov(X_c), C0, 5e-2);
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
        
        function verify_mp(g, d, m, cf, mu, C)
            assert(strcmp(g.tag, 'gaussd'));
            assert(isequal(g.ty, 'm'));
            assert(g.d == d);
            assert(g.n == m);
            assert(isequal(g.C.ty, cf));
            if ~isempty(mu)
                assert(isequal(g.mu, mu));
            end
            if ~isempty(C)
                assert(isequal(g.C, C));
            end
        end
        
        
        function verify_cp(g, d, m, cf, h, J)
            assert(strcmp(g.tag, 'gaussd'));
            assert(isequal(g.ty, 'c'));
            assert(g.d == d);
            assert(g.n == m);
            assert(isequal(g.J.ty, cf));
            if ~isempty(h)
                assert(isequal(g.h, h));
            end
            if ~isempty(J)
                assert(isequal(g.J, J));
            end
        end
        
        
        function v = calc_entropy(C)
            d = C.d;
            ldv = pdmat_lndet(C);
            v = 0.5 * (d * (log(2*pi) + 1) + ldv);
        end
        
        function v = calc_ca(mu, C)
            C = pdmat_fullform(C);
            v = mu' * (C \ mu);
        end
        
        function v = calc_cb(C)
            d = C.d;
            ldv = pdmat_lndet(C);
            v = -0.5 * (d * log(2*pi) + ldv);
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

        
        function G = gmle(X, w, cf, cov_tied)
            % A slow (but straightforward) way to implement MLE
            
            mu = vecmean(X, w);
            d = size(X, 1);
            m = size(w, 1);
            
            if cov_tied && m > 1
                sw = sum(w, 2);
                sw = sw / sum(sw);
            end
            
            switch cf
                case {'s', 'd'}
                    v = vecvar(X, w, mu);
                    if cf == 's'
                        v = mean(v, 1);
                    end
                    if cov_tied && m > 1
                        v = v * sw;
                    end
                    C = pdmat(cf, d, v);
                    
                case 'f'
                    v = veccov(X, w, mu);
                    if cov_tied && m > 1
                        v = reshape(v, d * d, m) * sw;
                        v = reshape(v, d, d);
                    end
                    C = pdmat(cf, d, v);
            end
            
            G = gaussd('m', mu, C);            
        end
        
    end
            
end
        

