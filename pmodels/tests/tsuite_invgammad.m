classdef tsuite_invgammad
    % Test suite for invgammad
    %
    
    % Create by Dahua Lin, on Sep 1, 2011
    %
    
    
    
    %% Properties
    
    properties
        dims = [1 2 5];
        nums = [1 3];
    end
    
    
    %% Test cases
    
    methods
        
        function test_basics(obj)
            run_multi(obj, @tsuite_invgammad.do_test_basics);
        end
        
        function test_statistics(obj)
            run_multi(obj, @tsuite_invgammad.do_test_statistics);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_invgammad.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_invgammad.do_test_sampling);
        end
    end
            
    
    %% Test Implementation
    
    methods(Static, Access='private')
        
        function do_test_basics(d, m, uniform_shape, share_scale)
            
            [g, A, B, gp] = tsuite_invgammad.make_obj(d, m, ...
                uniform_shape, share_scale);
            
            assert(isempty(g.lpconst));
            assert(isequal(size(gp.lpconst), [1, m]));
            
            lpc0 = bsxfun(@times, A, log(B)) - gammaln(A);
            lpc0 = sum(lpc0, 1);
            assert(isequal(size(lpc0), [1, m]));
            
            devcheck('lpconst calc', gp.lpconst, lpc0, 1e-12);
        end
    
        
        function do_test_statistics(d, m, uniform_shape, share_scale)
            
            [g, A, B] = tsuite_invgammad.make_obj(d, m, ...
                uniform_shape, share_scale);
            
            mean0 = bsxfun(@rdivide, B, A-1);
            var0 = bsxfun(@rdivide, B.^2, (A-1).^2 .* (A-2));
            mode0 = bsxfun(@rdivide, B, A+1);
            
            E1 = A + gammaln(A) - (1 + A) .* psi(A);
            E = bsxfun(@plus, E1, log(B));
            ent0 = sum(E, 1);
            
            mean1 = mean(g);
            var1 = var(g);
            mode1 = mode(g);
            ent1 = entropy(g);
            
            assert(isequal(size(mean1), [d, m]));
            assert(isequal(size(var1), [d, m]));
            assert(isequal(size(mode1), [d, m]));
            assert(isequal(size(ent1), [1, m]));
            
            devcheck('mean calc', mean0, mean1, 1e-14);
            devcheck('var calc',  var0, var1, 1e-14);
            devcheck('mode calc',  mode0, mode1, 1e-14);
            devcheck('entropy calc', ent0, ent1, 1e-12);
        end
    
        
        function do_test_evaluation(d, m, uniform_shape, share_scale)
            
            [g, A, B, gp] = tsuite_invgammad.make_obj(d, m, ...
                uniform_shape, share_scale);
            
            N = 100;
            X = rand(d, N) + 0.5;
            
            L0 = tsuite_invgammad.my_calc_logpdf(A, B, X);
            L1 = g.logpdf(X);
            L1p = gp.logpdf(X);
            assert(isequal(size(L1), [m, N]));
            assert(isequal(L1, L1p));
            
            L2 = zeros(m, N);
            for k = 1 : m
                L2(k, :) = gp.logpdf(X, k);
            end
            
            devcheck('logpdf eval', L1, L0, 1e-10);
            devcheck('logpdf eval (per-row)', L2, L1, 1e-13);
            
            P1 = g.pdf(X);
            assert(isequal(P1, exp(L1)));            
        end
        
        
        function do_test_sampling(d, m, uniform_shape, share_scale)
            
            g = tsuite_invgammad.make_obj(d, m, ...
                uniform_shape, share_scale);            
            
            mean0 = mean(g);
            var0 = var(g);
            
            ns = 1e6;
            if m == 1
                X1 = g.sample(ns);
                assert(isequal(size(X1), [d, ns]));
                
                devcheck('sample 1 - mean', vecmean(X1), mean0, 1e-2);
                devcheck('sample 1 - var',  vecvar(X1), var0, 8e-2);
            end
            
            X2 = g.sample(ns(ones(1, m)), 1:m);
            for k = 1 : m
                cX2 = X2(:, (k-1)*ns+1 : (k-1)*ns+ns);
                
                devcheck('sample 2 - mean', vecmean(cX2), mean0(:,k), 1e-2);
                devcheck('sample 2 - var',  vecvar(cX2), var0(:,k), 8e-2);
            end
        end
    end

    
    %% Auxiliary functions
    
    methods(Access='private')
        
        function run_multi(obj, tfunc)
            % run multiple test under different settings
            
            ds = obj.dims;
            ms = obj.nums;
            
            for d = ds
                for m = ms
                    tfunc(d, m, 0, 0);
                    tfunc(d, m, 0, 1);
                    tfunc(d, m, 1, 0);
                    tfunc(d, m, 1, 1);
                end
            end
            
        end
    end
    
    
    methods(Static, Access='private')
                
        function [g, A, B, gp] = make_obj(d, m, uniform_shape, share_scale)
            
            if uniform_shape
                alpha = rand(1, m) + 3;
                A = repmat(alpha, [d, 1]);
            else
                alpha = rand(d, m) + 3;
                A = alpha;
            end
            
            if share_scale
                beta = rand() + 0.5;
                B = repmat(beta, 1, m);
            else
                beta = rand(1, m) + 0.5;
                B = beta;
            end
            
            if ~uniform_shape
                g  = invgammad(alpha, beta);
                gp = invgammad(alpha, beta, [], 'pre');
            else
                g  = invgammad(alpha, beta, d);
                gp = invgammad(alpha, beta, d, 'pre');
            end
            
            assert(g.dim == d);
            assert(g.num == m);
            assert(isequal(g.alpha, alpha));
            assert(isequal(g.beta, beta));
            
            assert(gp.dim == d);
            assert(gp.num == m);
            assert(isequal(gp.alpha, alpha));
            assert(isequal(gp.beta, beta));
        end
                
        
        function L = my_calc_logpdf(A, B, X)
            
            [d, m] = size(A);
            n = size(X, 2);
            
            L = zeros(m, n);
            
            for k = 1 : m
                
                b = B(k);
                
                Pk = zeros(d, n);
                for i = 1 : d
                    a = A(i, k);
                    x = X(i, :);
                    Pk(i, :) = ((b./x).^a .* exp(-b./x)) ./ (x * gamma(a));
                end
                
                L(k, :) = sum(log(Pk), 1);
            end            
        end

    end
        
end
    



