classdef tsuite_dirichletd
    % Test suite for dirichlet distribution functions
    %
    
    %   History
    %   -------
    %       - Create by Dahua Lin, on Sep 1, 2011
    %       - Modified by Dahua Lin, on Sep 3, 2011
    %
    
    
    %% Properties
    
    properties        
        dims = [2 5];
        nums = [1 3];        
    end

    %% Test cases
    
    methods   
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_dirichletd.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_dirichletd.do_test_sampling);
        end
    end
    
    
    
    %% Test implementation
    
    methods(Static, Access='private')           
        
        function do_test_evaluation(d, m, alpha)
                        
            N = 100;
            X = rand(d, N);
            X = bsxfun(@times, X, 1 ./ sum(X, 1));
            
            if d > size(alpha, 1)
                A = repmat(alpha, d, 1);
            else
                A = alpha;
            end
            
            % covariance
            
            if m == 1
                C0 = tsuite_dirichletd.calc_cov(A);
                C1 = dird_cov(alpha, d);
                
                assert(isequal(size(C0), [d d]));
                assert(isequal(size(C1), [d d]));
                
                devcheck('cov', C1, C0, 1e-15);
            end
            
            % entropy
            
            ent0 = tsuite_dirichletd.calc_entropy(A);
            ent1 = dird_entropy(alpha, d);
            
            assert(isequal(size(ent0), [1 m]));
            assert(isequal(size(ent1), [1 m]));
            
            devcheck('entropy', ent1, ent0, 1e-12);            
            
            % log pdf
            
            L0 = tsuite_dirichletd.calc_logpdf(A, X);
            L1 = dird_logpdf(alpha, X);
            
            assert(isequal(size(L0), [m N]));
            assert(isequal(size(L1), [m N]));
            
            devcheck('logpdf', L1, L0, 1e-12);            
        end
        
        
        function do_test_sampling(d, m, alpha)
                        
            if m > 1
                return;
            end                
            
            if d > size(alpha, 1)
                A = repmat(alpha, d, 1);
            else
                A = alpha;
            end            
            
            mean0 = bsxfun(@times, A, 1 ./ sum(A, 1));
            cov0 = dird_cov(A);
            
            ns = 1e6;

            X = dird_sample(alpha, [d ns]);
            assert(isequal(size(X), [d, ns]));
            
            smean = vecmean(X);
            scov = veccov(X);
                        
            devcheck('sample - mean', smean, mean0, 1e-2);
            devcheck('sample - cov', scov, cov0, 2e-2);
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
                    
                    a0 = 1.2 + rand(1, m);
                    tfunc(d, m, a0);
                    
                    a1 = 1.2 + rand(d, m);
                    tfunc(d, m, a1);
                end
            end
            
        end
    end    
    
    
    methods(Static, Access='private')
        
        function C = calc_cov(A)
            
            d = size(A, 1);
            a0 = sum(A, 1);
            
            C = zeros(d, d);
            
            for i = 1 : d
                for j = 1 : d
                    
                    ai = A(i);
                    aj = A(j);
                    
                    if i == j
                        cv = ai * (a0 - ai) / (a0^2 * (a0 + 1));
                    else
                        cv = - ai * aj / (a0^2 * (a0 + 1));
                    end
                    
                    C(i, j) = cv;
                end
            end
        end
        
        
        function v = calc_entropy(A)
            
            logB = sum(gammaln(A), 1) - gammaln(sum(A, 1)); 
            
            d = size(A, 1);
            a0 = sum(A, 1);
            
            v = logB + (a0 - d) .* psi(a0) - sum((A - 1) .* psi(A), 1);
        end        
        
        function L = calc_logpdf(A, X)
            
            m = size(A, 2);
            n = size(X, 2);
            L = zeros(m, n);
                        
            logB = sum(gammaln(A), 1) - gammaln(sum(A, 1)); 
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
    
