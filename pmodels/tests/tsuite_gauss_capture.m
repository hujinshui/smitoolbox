classdef tsuite_gauss_capture
    % The test suite for various observation capturing functions for Gauss
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
        
        function test_gaussgm_capture(obj)
            run_multi(obj, @tsuite_gauss_capture.do_test_gaussgm_capture);
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
        
        function do_test_gaussgm_capture(cf, d, K, n)
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
            
            % prepare arguments
            
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
            end
                            
            % run the functions
            
            if ~use_A
                [dh0, dJ0] = tsuite_gauss_capture.gaussgm_capture_gt(X, w, Jx, []);            
                [dh, dJ] = gaussgm_capture(X, w, Jx);
                
                if cf == 's'
                    [dh2, dJ2] = gaussgm_capture(X, w, Jx.v);                                    
                    assert(isequal(dh, dh2));
                    assert(isequal(dJ, dJ2));
                end
            else                
                [dh0, dJ0] = tsuite_gauss_capture.gaussgm_capture_gt(X, w, Jx, A);
                [dh, dJ] = gaussgm_capture(X, w, Jx, A);
            end
            
            % verify results
           
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










