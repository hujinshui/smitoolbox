classdef tsuite_scale_invchi2d
    % Test suite for invchi2d
    %
    
    % Create by Dahua Lin, on Sep 1, 2011
    %

    
    %% Test cases
    
    methods
        
        function test_basics(obj)
            run_multi(obj, @tsuite_scale_invchi2d.do_test_basics);
        end
        
        function test_statistics(obj)
            run_multi(obj, @tsuite_scale_invchi2d.do_test_statistics);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_scale_invchi2d.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_scale_invchi2d.do_test_sampling);
        end
    end
    
    
    %% Test implementation
    
    methods(Static, Access='private')
        
        function do_test_basics(d, is_uniform)
            
            D = tsuite_scale_invchi2d.make_obj(d, is_uniform);
            
            alpha0 = D.nu / 2;
            beta0 = D.nu * D.sigma2 / 2;
            
            assert(norm(D.alpha - alpha0) < 2e-15);
            assert(norm(D.beta - beta0) < 2e-15);
        end
        
        
        function do_test_statistics(d, is_uniform)
            
            D = tsuite_scale_invchi2d.make_obj(d, is_uniform);
            nu = D.nu;
            sigma2 = D.sigma2;
            
            mean0 = (nu * sigma2) / (nu - 2);
            var0 = (2 * nu^2 .* (sigma2.^2)) / ((nu - 2)^2 * (nu - 4));
            mode0 = (nu * sigma2) / (nu + 2);
            
            if is_uniform && d > 1
                mean0 = mean0(ones(d,1), :);
                var0 = var0(ones(d,1), :);
                mode0 = mode0(ones(d,1), :);
            end
            
            if is_uniform
                ent0 = d * (invgamma_entropy(D.alpha) + log(D.beta));
            else
                ent0 = d * invgamma_entropy(D.alpha) + sum(log(D.beta));
            end
            
            assert(isequal(size(mean(D)), [d 1]));
            assert(isequal(size(var(D)), [d 1]));
            assert(isequal(size(mode(D)), [d 1]));
            assert(isscalar(entropy(D)));
            
            devcheck('mean calc', mean(D), mean0, 1e-13);
            devcheck('var calc', var(D), var0, 1e-13);
            devcheck('mode calc', mode(D), mode0, 1e-13);
            devcheck('entropy', entropy(D), ent0, 1e-13);
        end
        
        
        function do_test_evaluation(d, is_uniform)
            
            [D, Sig2] = tsuite_scale_invchi2d.make_obj(d, is_uniform);
            
            n = 100;
            X = rand(d, n) + 1;
            
            L0 = tsuite_scale_invchi2d.my_calc_logpdf(D.nu, Sig2, X);
            L = D.logpdf(X);
            assert(isequal(size(L), [1, n]));
            devcheck('logpdf', L, L0, 1e-12);
            
            P = D.pdf(X);
            assert(isequal(size(P), [1, n]));
            devcheck('pdf', P, exp(L), 1e-13);            
        end
        
        
        function do_test_sampling(d, is_uniform)
            
            D = tsuite_scale_invchi2d.make_obj(d, is_uniform);
            
            ns = 5e6;
            X = D.sample(ns);
            
            devcheck('sample mean', vecmean(X), mean(D), 2e-2);
            devcheck('sample var', vecvar(X), var(D), 5e-2);            
        end
    end
        
    
    %% Auxiliary functions
    
    methods(Access='private')
        
        function run_multi(obj, tfunc) %#ok<MANU>
            
            tfunc(1, false);
            tfunc(1, true);
            tfunc(3, false);
            tfunc(3, true);            
        end        
    end
    
    
    methods(Static, Access='private')    
        
        function [D, Sig2] = make_obj(d, is_uniform)
            
            nu = rand() + 8;
            
            if is_uniform
                sigma2 = rand() + 2;
            else
                sigma2 = rand(d, 1) + 2;
            end
            
            if is_uniform
                D = scale_invchi2d(nu, sigma2, d);
                Sig2 = sigma2(ones(d, 1), 1);
            else
                D = scale_invchi2d(nu, sigma2);
                Sig2 = sigma2;
            end
            
            assert(D.num == 1);
            assert(D.dim == d);
            assert(isequal(D.nu, nu));
            assert(isequal(D.sigma2, sigma2));            
        end
        
        function L = my_calc_logpdf(nu, sigma2, X)
            
            d = size(sigma2, 1);
            L = zeros(d, size(X,2));
            for i = 1 : d
                x = X(i, :);
                c = (nu / 2) * log(sigma2(i) * nu / 2) - gammaln(nu / 2);
                t1 = (nu * sigma2(i)) ./ (2 * x);
                t2 = (1 + nu / 2) * log(x);
                L(i, :) = c - t1 - t2;
            end
            L = sum(L, 1);
        end
    
    end
    
end
    
