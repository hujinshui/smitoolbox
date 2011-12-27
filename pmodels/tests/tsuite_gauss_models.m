classdef tsuite_gauss_models
    % The test suite for various Gauss based generative models
    %
    
    % Created by Dahua Lin, on Dec 13, 2011
    %
    
    %% configurations
    
    properties
        Jx_forms = {'s', 'd', 'f'};
        dims = {1, 2, 5, [3, 2]};
        Ks = [1, 3];
    end
    
    %% main test cases
    
    methods
        
        function test_gaussgen(obj)
            run_multi(obj, @tsuite_gauss_models.do_test_gaussgen);
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
        
        function do_test_gaussgen(cf, d, K, n)
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
                        
            X = randn(d, n);
            U = randn(q, K);
            
            if use_A
                A = randn(d, q);
            else
                A = [];
            end
            
            if ~use_A
                g0 = gaussgen(d);
                gm = gaussgen(Jx);
            else
                g0 = gaussgen(d, q);
                gm = gaussgen(Jx, A);
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
            
            assert(g0.query_params(U) == K);
            assert(gm.query_params(U) == K);
            
            % verify loglik evaluation
            
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
            
            devcheck('simplegen (LL)', LL, LL0, 1e-12);            
            
            % verify capturing (conjugate updates)
            
            if K == 1
                Zi = [];
            else
                ssiz = floor(n / K);
                Zi = cell(1, K);
                for k = 1 : K
                    Zi{k} = (k-1)*ssiz + 1 : k * ssiz;
                end
            end
            Zw = rand(K, n);                        
            
            [dh0_i, dJ0_i] = tsuite_gauss_models.gaussgen_capture_gt(X, Zi, Jx, A);
            dG_i = gm.capture(X, Zi);
            
            [dh0_w, dJ0_w] = tsuite_gauss_models.gaussgen_capture_gt(X, Zw, Jx, A);
            dG_w = gm.capture(X, Zw);
           
            assert(isequal(size(dh0_i), [q, K]));
            assert(is_pdmat(dJ0_i) && dJ0_i.d == q && dJ0_i.n == K);            
            
            assert(is_gaussd(dG_i) && dG_i.ty == 'c' && dG_i.d == q && dG_i.n == K);
            assert(isequal(size(dG_i.h), [q, K]));
            assert(dG_i.J.d == q && dG_i.J.n == K);
            
            assert(is_gaussd(dG_w) && dG_w.ty == 'c' && dG_w.d == q && dG_w.n == K);
            assert(isequal(size(dG_w.h), [q, K]));
            assert(dG_w.J.d == q && dG_w.J.n == K);
                        
            devcheck('gaussgen w/ Zi (dh)', dG_i.h, dh0_i, 1e-12);
            devcheck('gaussgen w/ Zi (dJ)', dG_i.J.v, dJ0_i.v, 1e-12);
            
            devcheck('gaussgen w/ Zw (dh)', dG_w.h, dh0_w, 1e-12);
            devcheck('gaussgen w/ Zw (dJ)', dG_w.J.v, dJ0_w.v, 1e-12);
            
            
            % verify MLE
            
            Ui = gm.mle(X, Zi);
            Ui0 = pdmat_lsolve(dJ0_i, dh0_i);
            
            Uw = gm.mle(X, Zw);
            Uw0 = pdmat_lsolve(dJ0_w, dh0_w);
            
            assert(isequal(size(Ui0), [q, K]));
            assert(isequal(size(Ui), [q, K]));
            devcheck('gaussgen (Ui)', Ui, Ui0, 1e-10);  
                                    
            assert(isequal(size(Uw0), [q, K]));
            assert(isequal(size(Uw), [q, K]));
            devcheck('gaussgen (Uw)', Uw, Uw0, 1e-10);  
        end
        
        
        function [dh, dJ] = gaussgen_capture_gt(X, Z, Jx, A)
            % Calculate the ground-truth for gaussgm_capture
            
            n = size(X, 2);            
            [zty, K] = verify_Zarg(Z, n);
            
            if zty == 0
                w = ones(1, n);
            elseif zty == 1
                w = Z;
            elseif zty == 2
                w = zeros(K, n);
                for k = 1 : K
                    w(k, Z{k}) = 1;
                end
            end
            
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










