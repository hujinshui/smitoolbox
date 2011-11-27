classdef tsuite_dirichletd
    % Test suite for gammd class
    %
    
    %   History
    %   -------
    %       - Create by Dahua Lin, on Sep 1, 2011
    %       - Modified by Dahua Lin, on Sep 3, 2011
    %
    
    
    %% Properties
    
    properties        
        Ks = [1 2 5];
        nums = [1 3];        
    end

    %% Test cases
    
    methods
    
        function test_basics(obj)
            run_multi(obj, @tsuite_dirichletd.do_test_basics);
        end
        
        function test_statistics(obj)
            run_multi(obj, @tsuite_dirichletd.do_test_statistics);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_dirichletd.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_dirichletd.do_test_sampling);
        end
    end
    
    
    
    %% Test implementation
    
    methods(Static, Access='private')
    
        function do_test_basics(K, m, is_sym)
            
            [D, Alpha, Dp] = tsuite_dirichletd.make_obj(K, m, is_sym);  %#ok<ASGLU>

            logB0 = tsuite_dirichletd.my_calc_logB(Alpha);
            assert(max(abs(logB0 - Dp.logB)) < 1e-14);
        end
        
        
        function do_test_statistics(K, m, is_sym)
            
            [D, A, Dp] = tsuite_dirichletd.make_obj(K, m, is_sym);
            
            % compute ground truths
            
            sumA = sum(A, 1);
            
            mean0 = bsxfun(@times, A, 1 ./ sumA);
            
            var0 = zeros(K, m);
            for i = 1 : m
                for k = 1 : K
                    a = A(k, i);
                    a0 = sumA(i);
                    var0(k, i) = a * (a0 - a) / (a0^2 * (a0 + 1));
                end
            end
            
            if m == 1
                cov0 = zeros(K, K);
                a0 = sumA;
                for k = 1 : K
                    for l = 1 : K
                        if k == l
                            cv = A(k) * (a0 - A(k)) / (a0^2 * (a0 + 1));
                        else
                            cv = - A(k) * A(l) / (a0^2 * (a0 + 1));
                        end
                        cov0(k, l) = cv;
                    end                    
                end
            end
                                                
            mode0 = bsxfun(@times, A-1, 1 ./ (sumA - K));
            
            logB = Dp.logB;            
            ent0 = bsxfun(@minus, logB + (sumA - K) .* psi(sumA), dot(A-1, psi(A), 1));
            
            % compare with the results
                        
            mean1 = mean(D);
            var1 = var(D);
            if m == 1
                cov1 = cov(D);
            end            
            mode1 = mode(D);
            ent1 = entropy(D);
                                    
            mean2 = mean(Dp);
            var2 = var(Dp);
            if m == 1
                cov2 = cov(Dp);
            end            
            mode2 = mode(Dp);
            ent2 = entropy(Dp);
            
            assert(isequal(size(mean1), [K, m]));
            assert(isequal(size(var1), [K, m]));
            if m == 1
                assert(isequal(size(cov1), [K, K]));
            end
            assert(isequal(size(mode1), [K, m]));
            assert(isequal(size(ent1), [1, m]));
            
            assert(isequal(size(mean2), [K, m]));
            assert(isequal(size(var2), [K, m]));
            if m == 1
                assert(isequal(size(cov2), [K, K]));
            end
            assert(isequal(size(mode2), [K, m]));
            assert(isequal(size(ent2), [1, m]));
                        
            devcheck('mean1 calc', mean0, mean1, 1e-14);
            devcheck('var1 calc',  var0, var1, 1e-14);
            if m == 1
                devcheck('cov1 calc', cov0, cov1, 1e-14);
            end
            devcheck('mode1 calc',  mode0, mode1, 1e-14);
            devcheck('entropy1 calc', ent0, ent1, 1e-12);
            
            devcheck('mean2 calc', mean0, mean2, 1e-14);
            devcheck('var2 calc',  var0, var2, 1e-14);
            if m == 1
                devcheck('cov2 calc', cov0, cov2, 1e-14);
            end
            devcheck('mode2 calc',  mode0, mode2, 1e-14);
            devcheck('entropy2 calc', ent0, ent2, 1e-12);
        end
        
        
        function do_test_evaluation(K, m, is_sym)
            
            [D, A, Dp] = tsuite_dirichletd.make_obj(K, m, is_sym);
            
            N = 100;
            X = rand(K, N);
            X = bsxfun(@times, X, 1 ./ sum(X, 1));
            
            L0 = tsuite_dirichletd.my_calc_logpdf(A, X);
            L1 = D.logpdf(X);
            L1p = Dp.logpdf(X);
            assert(isequal(size(L1), [m, N]));
            assert(isequal(L1, L1p));
            
            L2 = zeros(m, N);
            for k = 1 : m
                L2(k, :) = Dp.logpdf(X, k);
            end
            
            devcheck('logpdf eval', L1, L0, 1e-12);
            devcheck('logpdf eval (per-row)', L2, L1, 1e-13);
            
            P1 = D.pdf(X);
            P1p = Dp.pdf(X);
            assert(isequal(P1, exp(L1)));            
            assert(isequal(P1p, exp(L1p)));
        end
        
        
        function do_test_sampling(K, m, is_sym)
            
            g = tsuite_dirichletd.make_obj(K, m, is_sym);
            
            mean0 = mean(g);
            var0 = var(g);
            
            ns = 5e5;
            if m == 1
                X1 = g.sample(ns);
                assert(isequal(size(X1), [K, ns]));
                
                devcheck('sample 1 - mean', vecmean(X1), mean0, 2e-2);
                devcheck('sample 1 - var',  vecvar(X1), var0, 5e-2);
            end
            
            X2 = g.sample(ns(ones(1, m)), 1:m);
            assert(isequal(size(X2), [K, ns * m]));
            for k = 1 : m
                cX2 = X2(:, (k-1)*ns+1 : (k-1)*ns+ns);
                
                devcheck('sample 2 - mean', vecmean(cX2), mean0(:,k), 2e-2);
                devcheck('sample 2 - var',  vecvar(cX2), var0(:,k), 5e-2);
            end
        end
    end
    
    
    
    %% Auxiliary functions
    
    methods(Access='private')
        
        function run_multi(obj, tfunc)
            % run multiple test under different settings
            
            ks = obj.Ks;
            ms = obj.nums;
            
            for k = ks
                for m = ms
                    tfunc(k, m, 0);
                    tfunc(k, m, 1);
                end
            end
            
        end
    end    
    
    
    methods(Static, Access='private')
        
        function [D, Alpha, Dp] = make_obj(K, m, is_sym)
        
            if is_sym
                alpha = 1.2 + rand(K, m);
                Alpha = alpha;
            else
                alpha = 1.2 + rand(1, m);
                Alpha = repmat(alpha, K, 1);
            end
            
            D = dirichletd(K, alpha);
            Dp = dirichletd(K, alpha, 'pre');            
            
            assert(D.K == K);
            assert(D.num == m);
            assert(isequal(D.alpha, alpha));
            assert(isempty(D.logB));
            
            assert(Dp.K == K);
            assert(Dp.num == m);
            assert(isequal(Dp.alpha, alpha));
            assert(isequal(size(Dp.logB), [1 m]));
        end
        
        
        function v = my_calc_logB(Alpha)
            v = sum(gammaln(Alpha), 1) - gammaln(sum(Alpha, 1));             
        end
        
        function L = my_calc_logpdf(A, X)
            
            m = size(A, 2);
            n = size(X, 2);
            L = zeros(m, n);
                        
            logB = tsuite_dirichletd.my_calc_logB(A);
            for i = 1 : m
                a = A(:, i);
                for j = 1 : n
                    x = X(:, j);
                    
                    L(i, j) = -logB(i) + (a-1)' * log(x);
                end
            end
        end
        
    end    
        
end
    
