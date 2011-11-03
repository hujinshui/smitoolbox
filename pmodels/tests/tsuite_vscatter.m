classdef tsuite_vscatter
    % A test suite for vscatter function
    %
    
    % Created by Dahua Lin, on Sep 28, 2011
    %
    
    %% Test cases
    
    methods
        
        function test_scatter_vec(obj) %#ok<MANU>
            tsuite_vscatter.run_multi(@tsuite_vscatter.do_test_scatter_vec);
        end
        
        function test_scatter_mat(obj) %#ok<MANU>
            tsuite_vscatter.run_multi(@tsuite_vscatter.do_test_scatter_mat);
        end    
    end
    
    
    methods(Static)
        
        function do_test_scatter_vec(X, U, W)
            
            [d, n] = size(X);
            K = size(U, 2);
            
            R0 = zeros(d, K);
            for k = 1 : K
                if isequal(U, 0)
                    u = zeros(d, 1);
                else
                    u = U(:,k);
                end
                if isscalar(W)
                    w = W * ones(1, n);
                else
                    w = W(k,:);
                end
                
                Y = bsxfun(@minus, X, u);
                R0(:,k) = (Y.^2) * w';
            end
            
            R = vscatter(X, U, W, 'v');
            assert(isequal(size(R0), size(R)));
            assert(max(abs(R(:) - R0(:))) < 1e-12);
        end
        
        
        function do_test_scatter_mat(X, U, W)
            
            [d, n] = size(X);
            K = size(U, 2);
            
            R0 = zeros(d, d, K);
            for k = 1 : K
                if isequal(U, 0)
                    u = zeros(d, 1);
                else
                    u = U(:,k);
                end
                if isscalar(W)
                    w = W * ones(1, n);
                else
                    w = W(k,:);
                end
                
                Y = bsxfun(@minus, X, u);
                C = Y * diag(w) * Y';
                
                R0(:,:,k) = 0.5 * (C + C');
            end
            
            R = vscatter(X, U, W, 'c');                        
            assert(isequal(size(R0), size(R)));            
            for k = 1 : K
                C = R(:,:,k);
                assert(isequal(C, C'));
            end
            assert(max(abs(R(:) - R0(:))) < 1e-12);
        end
        
        
        function run_multi(tfunc)
            
            n = 50;
            
            for d = [1 2 5]
                for K = [1 3]
                    X = rand(d, n);
                    U = rand(d, K);
                    W = rand(K, n);
                    
                    if K == 1
                        tfunc(X, 0, 1);
                        tfunc(X, 0, 2.5);
                        tfunc(X, 0, W);
                    end
                    
                    tfunc(X, U, 1);
                    tfunc(X, U, 2.5);
                    tfunc(X, U, W);                    
                end
            end
        end
        
    end
    
    
end
