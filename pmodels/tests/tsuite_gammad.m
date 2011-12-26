classdef tsuite_gammad
    % Test suite for gammd class
    %
    
    %   History
    %   -------
    %       - Create by Dahua Lin, on Sep 1, 2011
    %       - Modified by Dahua Lin, on Sep 3, 2011
    %
    
    
    %% Properties
    
    properties        
        dims = { [1 1 1], [1 1 3], [1 3 3], [3 1 3], [3 3 3] };
        nums = { [1 1], [1 5], [5 1], [5 5] };        
    end

    %% Test cases
    
    methods
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_gammad.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_gammad.do_test_sampling);
        end
    end
    
    
    
    %% Test implementation
    
    methods(Static, Access='private')    
                
        function do_test_evaluation(d, m, alpha, beta)
                        
            N = 100;
            X = rand(d, N);
            
            % entropy
            
            ent0 = tsuite_gammad.calc_entropy(alpha, beta);
            ent = gammad_entropy(alpha, beta);
            
            assert(isequal(size(ent0), [1 m]));
            assert(isequal(size(ent), [1 m]));
            devcheck('entropy eval', ent, ent0, 1e-12);
            
            % logpdf
            
            L0 = tsuite_gammad.calc_logpdf(alpha, beta, X);
            L = gammad_logpdf(alpha, beta, X);
            
            assert(isequal(size(L0), [m N]));
            assert(isequal(size(L), [m N]));            
                        
            devcheck('logpdf eval', L, L0, 1e-12);          
        end
        
        
        function do_test_sampling(d, m, alpha, beta)
            
            if m > 1
                return;
            end                        
            
            mean0 = alpha .* beta;
            var0 = alpha .* (beta.^2);
            
            if d > 1 && isscalar(alpha) && isscalar(beta)
                mean0 = mean0(ones(d, 1), 1);
                var0 = var0(ones(d, 1), 1);
            end
            
            ns = 1e5;
            X1 = gammad_sample(alpha, beta, [d ns]);
            assert(isequal(size(X1), [d, ns]));
            
            devcheck('sample 1 - mean', vecmean(X1), mean0, 2e-2);
            devcheck('sample 1 - var',  vecvar(X1), var0, 0.15);
        end
    end
    
    
    
    %% Auxiliary functions
    
    methods(Access='private')
        
        function run_multi(obj, tfunc)
            % run multiple test under different settings
            
            ds = obj.dims;
            ms = obj.nums;
                        
            for i = 1 : numel(ds)
                for j = 1 : numel(ms)
                    
                    d = ds{i};
                    m = ms{j};
                    
                    da = d(1);
                    db = d(2);
                    d = d(3);
                    
                    ma = m(1);
                    mb = m(2);
                    m = max(ma, mb);
                    
                    alpha = rand(da, ma) + 1.5;
                    beta = rand(db, mb) + 0.5;
                    
                    tfunc(d, m, alpha, beta);
                end
            end
            
        end
    end    
    
    
    methods(Static, Access='private')                
        
        function v = calc_entropy(A, B)            
            v = bsxfun(@plus, A + gammaln(A) + (1 - A) .* psi(A), log(B));
            v = sum(v, 1);
        end
        
        function L = calc_logpdf(A, B, X)
            
            m = max(size(A, 2), size(B, 2));
            [d, n] = size(X);
            
            A = bsxfun(@times, ones(d, m), A);
            B = bsxfun(@times, ones(d, m), B);
                      
            L = zeros(m, n);
            
            for k = 1 : m                
                Pk = zeros(d, n);
                for i = 1 : d
                    a = A(i, k);
                    b = B(i, k);
                    Pk(i, :) = gampdf(X(i, :), a, b);
                end
                
                L(k, :) = sum(log(Pk), 1);
            end
        end
        
    end    
        
end
    
