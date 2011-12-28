classdef tsuite_ppca
    % A test suite for PPCA functions
    %
    
    % Created by Dahua Lin, on Dec 27, 2011
    %
    
    %% properties
    
    properties
        dims = {[2 1], [5, 3], [10, 4]};                
    end
    
    %% Test cases
    
    methods    
        function test_basics(obj)
            run_multi(obj, @tsuite_ppca.do_test_basics);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_ppca.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_ppca.do_test_sampling);
        end        
    end
    
    %% Core Test functions
    
    methods(Static)
        
        function do_test_basics(d, q, zmean)
            % Test basic functionalities
            
            M = tsuite_ppca.rand_model(d, q, zmean);
            
            % check W
            
            W0 = tsuite_ppca.Wmat(M);            
            W = ppca_W(M);
            
            assert(isequal(size(W), [d q]));
            devcheck('Wmat', W, W0, 1e-15);
            
            % check cov and icov
            
            C0 = tsuite_ppca.Cmat(M);
            C = ppca_cov(M);
            
            assert(isa(C, 'double') && isequal(size(C), [d d]));
            devcheck('Cov', C, C0, 1e-13);
            
            J = ppca_icov(M);
            assert(isa(J, 'double') && isequal(size(J), [d d]));
            
            devcheck('ICov', C * J, eye(d), 1e-12);
            devcheck('ICov', J * C, eye(d), 1e-12);
            
            % check conversion to gaussd
            
            Gm = gaussd('m', M);
            Gc = gaussd('c', M);
            
            assert(is_gaussd(Gm) && Gm.d == d && Gm.n == 1);
            assert(is_gaussd(Gc) && Gc.d == d && Gc.n == 1);
            
            assert(is_pdmat(Gm.C) && Gm.C.d == d && Gm.C.ty == 'f' && Gm.C.n == 1);
            assert(is_pdmat(Gc.J) && Gc.J.d == d && Gc.J.ty == 'f' && Gc.J.n == 1);
            
            assert(isequal(Gm.C.v, C));
            assert(isequal(Gc.J.v, J));
            
            if zmean
                assert(isequal(Gm.mu, 0));
                assert(isequal(Gc.h, 0));
            else
                assert(isequal(size(Gm.mu), [d 1]));
                assert(isequal(size(Gc.h), [d 1]));
                
                assert(isequal(Gm.mu, M.mu));
                
                h0 = Gm.C.v \ Gm.mu;
                devcheck('h-vec', Gc.h, h0, 1e-13);
            end                                
        end
        
        
        function do_test_evaluation(d, q, zmean)
            % Test evaluation and transform
            
            M = tsuite_ppca.rand_model(d, q, zmean);
            
            % check ft and bt
            
            n = 100;
            Z = randn(q, n);
            X0 = tsuite_ppca.ftr(M, Z);
            X = ppca_ft(M, Z);
            
            assert(isequal(size(X), [d, n]));
            devcheck('ft', X, X0, 1e-14);
            
            Zr = ppca_bt(M, X);
            
            assert(isequal(size(Z), [q n]));
            devcheck('bt', Zr, Z, 1e-12);
            
            % sqmahdist
            
            Gm = gaussd('m', M);
            Y = randn(d, n);
            
            D0 = gaussd_sqmahdist(Gm, Y);
            D = ppca_sqmahdist(M, Y);
            
            assert(isequal(size(D), [1 n]));
            devcheck('sqmahdist', D, D0, 1e-11);
            
            % logpdf
            
            L0 = gaussd_logpdf(Gm, Y);
            L = ppca_logpdf(M, Y);
            
            assert(isequal(size(L), [1 n]));
            devcheck('logpdf', L, L0, 1e-11);            
        end
              
        
        function do_test_sampling(d, q, zmean)
            % Test sampling
            
            M = tsuite_ppca.rand_model(d, q, zmean);
            
            C = ppca_cov(M);
            
            if zmean
                mu = zeros(d, 1);
            else
                mu = M.mu;
            end
            
            ns = 1e6;
            X = ppca_sample(M, ns);
            
            Xmean = vecmean(X);
            Xcov = veccov(X);
            
            devcheck('sample mean', Xmean, mu, 5e-2);
            devcheck('sample cov', Xcov, C, 2e-1);
        end
        
    end
    
    
    
    %% Auxiliary functions
    
    methods(Access='private')
        
        function run_multi(obj, tfunc)
            % Run a test case under multiple conditions
            
            ds = obj.dims;
            
            for i = 1 : numel(ds)
                d = ds{i}(1);
                q = ds{i}(2);
 
                tfunc(d, q, true);
                tfunc(d, q, false);
            end                        
        end
    end
    
    
    methods(Static, Access='private')
        
        function M = rand_model(d, q, zmean)
            
            U = orth(randn(d, d));
            U = U(:, 1:q);
            
            s = sort(rand(1, q) * 4 + 1);
            se = 0.5;
            
            if zmean
                mu = 0;
            else
                mu = randn(d, 1);
            end
            
            M = ppca_model(U, s, se, mu);
            
            assert(is_ppca(M));
            assert(M.d == d);
            assert(M.q == q);
            assert(isequal(M.s, s));
            assert(M.se == se);
            assert(isequal(M.B, U));
            assert(isequal(M.mu, mu));
        end
        
        
        function W = Wmat(M)
            W = bsxfun(@times, M.B, M.s);            
        end
        
        function C = Cmat(M)
            D = diag(M.s.^2);
            C = M.B * D * M.B';
            C = adddiag(C, M.se^2);
            if ~isequal(C, C')
                C = 0.5 * (C + C');
            end
        end
        
        
        function X = ftr(M, Z)            
            W = ppca_W(M);
            X = bsxfun(@plus, W * Z, M.mu);             
        end
        
    end
    
    
end

