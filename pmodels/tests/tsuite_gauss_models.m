classdef tsuite_gauss_models
    % The test suite for various Gauss based generative models
    %
    
    % Created by Dahua Lin, on Dec 13, 2011
    %
    
    %% configurations
    
    properties
        Jx_forms = {'s', 'd', 'f'};
        dims = {1, 2, 5, [3, 2]};
        Ks = [0, 1, 3];
    end
    
    %% main test cases
    
    methods
        
        function test_gaussgm(obj)
            run_multi(obj, @tsuite_gauss_models.do_test_gaussgm);
        end
                
        
        function run_multi(obj, tfunc)
            % Run a test case under multiple conditions
            
            cfs = obj.Jx_forms;
            ds = obj.dims;
            ks = obj.Ks;
            n = 25;
            
            for i = 1 : numel(cfs)
                cf = cfs{i};
                for j = 1 : numel(ds)
                    d = ds{j}; 
                    for k = ks
                        tfunc(cf, d, k, n);
                    end
                end
            end
            
        end
        
    end
                  
    
    %% core testing functions
    
    methods(Static, Access='private')
        
        %% gaussgn_capture
        
        function do_test_gaussgm(cf, d, K, n)
            % Perform the test of gaussgm_capture on a specific setting
            %
            %   cf:     the form of Jx
            %   d:      the dimensions
            %   K:      the number of rows in w
            %   n:      the number of observed samples
            %
            
            % parse inputs
            
            if isscalar(d)
                use_A = false;
                q = d;
            else
                use_A = true;
                q = d(2);
                d = d(1);
            end
            
            % prepare model
            
            Jx = rand_pdmat(cf, d, 1, [1 2]);
            
            if K == 0
                w = [];
                K = 1;
            else
                w = rand(K, n);
            end
            
            X = randn(d, n);
            
            if use_A
                A = randn(d, q);
            else
                A = [];
            end
            
            if ~use_A
                g0 = gaussgm(d);
                gm = gaussgm(Jx);
            else
                g0 = gaussgm(d, q);
                gm = gaussgm(Jx, A);
            end
            
            % verify models
            
            assert(g0.xdim == d);
            assert(g0.pdim == q);
            assert(isempty(g0.Gx));
            assert(isempty(g0.Gx_cb));
            assert(isempty(g0.Jx));
            assert(isempty(g0.A));
            
            assert(gm.xdim == d);
            assert(gm.pdim == q);
            assert(is_pdmat(gm.Jx) && isequal(gm.Jx, Jx));
            assert(isempty(gm.A) == ~use_A);
            if use_A
                assert(isequal(gm.A, A));
            end
            assert(is_gaussd(gm.Gx) && gm.Gx.ty == 'c' && gm.Gx.n == 1);
            assert(isequal(gm.Gx.h, 0));
            assert(isequal(gm.Gx.J, Jx));
            assert(isscalar(gm.Gx_cb));
            
            [~, cb] = gaussd_const(gm.Gx);
            devcheck('gaussgm (cb)', cb, gm.Gx_cb, 1e-15);                        
                            
            assert(g0.query_obs(X) == n);
            assert(gm.query_obs(X) == n);
            
            % verify loglik evaluation
            
            U = randn(q, K);
            
            LL = gm.loglik(U, X);
            
            if ~use_A
                AU = U;
            else
                AU = A * U;
            end            
            
            LL0 = zeros(K, n);
            for k = 1 : K
                lv = gaussd_logpdf(gm.Gx, bsxfun(@minus, X, AU(:,k)));
                LL0(k,:) = lv;
            end
            
            devcheck('gaussgm (LL)', LL, LL0, 1e-12);            
            
            % verify capturing (conjugate updates)
            
            [dh0, dJ0] = tsuite_gauss_models.gaussgm_capture_gt(X, w, Jx, A);
            [dh, dJ] = gm.capture(X, w);            
           
            assert(isequal(size(dh0), [q, K]));
            assert(isequal(size(dh), [q, K]));
            assert(is_pdmat(dJ0) && dJ0.d == q && dJ0.n == K);            
            assert(is_pdmat(dJ) && dJ.d == q && dJ.n == K);
            assert(isequal(dJ0.ty, dJ.ty));
            
            devcheck('gaussgm_capture (dh)', dh, dh0, 1e-12);
            devcheck('gaussgm_capture (dJ)', dJ.v, dJ0.v, 1e-12);
            
        end
        
        
        function [dh, dJ] = gaussgm_capture_gt(X, w, Jx, A)
            % Calculate the ground-truth for gaussgm_capture
            
            n = size(X, 2);
            if isempty(w)
                w = ones(1, n);
            end
            K = size(w, 1);
            
            if isempty(A)
                q = size(X, 1);
            else
                q = size(A, 2);
            end
            
            dh = zeros(q, K);
            for k = 1 : K
                dh_k = pdmat_mvmul(Jx, (X * w(k,:)'));
                if ~isempty(A)
                    dh_k = A' * dh_k;
                end
                dh(:, k) = dh_k;
            end
            
            sw = sum(w, 2).';
            if isempty(A)                
                dJ = pdmat_scale(Jx, sw);
            else
                dJ = zeros(q, q, K);
                Jx = pdmat_fullform(Jx);
                for k = 1 : K                    
                    dJ_k = sw(k) * (A' * Jx * A);
                    dJ_k = 0.5 * (dJ_k + dJ_k');
                    dJ(:,:,k) = dJ_k;
                end
                dJ = pdmat('f', q, dJ);
            end
            
            
        end               
        
        
    end
            
    
end










