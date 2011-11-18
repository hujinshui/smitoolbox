classdef tsuite_wishartd
    % Test suite for wishartd and wishart_sample
    %
    
    % Created by Dahua Lin, on Sep 3, 2011
    %
   
    
    %% Properties
    
    properties
        dims = [1 2 3 5];
        types = {'s', 'd', 'f'};
        deg = 8.5;        
    end
        
    %% Test cases
    
    methods
    
        function test_basics(obj)
            run_multi(obj, @tsuite_wishartd.do_test_basics);
        end
        
        function test_statistics(obj)
            run_multi(obj, @tsuite_wishartd.do_test_statistics);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_wishartd.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_wishartd.do_test_sampling);
        end
    end
    
    
    %% Test implementation
    
    methods(Static, Access='private')
        
        function do_test_basics(d, ty, deg)
            
            S = rand_pdmat(ty, d, 1, [0.5, 1.5]);
            
            D0 = wishartd(S, deg);
            D1 = wishartd(S, deg, 'pre');
            
            assert(D0.num == 1);
            assert(D0.dim == d);
            assert(D0.deg == deg);
            assert(isequal(D0.S, S));
            assert(isempty(D0.invS));
            assert(isempty(D0.lpconst));
            
            assert(D1.num == 1);
            assert(D1.dim == d);
            assert(D1.deg == deg);
            assert(isequal(D1.S, S));
            assert(~isempty(D1.invS));
            assert(~isempty(D1.lpconst));
                        
            q = deg / 2;
            lpc = q*d*log(2) + q * pdmat_lndet(S) + mvgammaln(d, q);
            
            devcheck('lpconst', D1.lpconst, lpc, 1e-13);
            assert(isequal(pdmat_inv(S), D1.invS));            
        end
        
        
        function do_test_statistics(d, ty, deg)
            
            S = rand_pdmat(ty, d, 1, [0.5, 1.5]);  
            Sm = pdmat_fullform(S);
            D0 = wishartd(S, deg);
            
            Mean0 = deg * Sm;
            Mode0 = (deg - d - 1) * Sm;
            
            devcheck('mean', pdmat_fullform(mean(D0)), Mean0, 1e-13);
            devcheck('mode', pdmat_fullform(mode(D0)), Mode0, 1e-13);
        end
    
        
        function do_test_evaluation(d, ty, deg)
            
            S = rand_pdmat(ty, d, 1, [0.5, 1.5]);
            Sm = pdmat_fullform(S);
            D0 = wishartd(S, deg);
            D1 = wishartd(S, deg, 'pre');
                        
            n = 100;
            W_s = rand_pdmat('s', d, n, [0.5 1.5]);
            W_d = rand_pdmat('d', d, n, [0.5 1.5]);
            W_f = rand_pdmat('f', d, n, [0.5 1.5]);
            
            L1_s0 = D0.logpdf(W_s);
            L1_d0 = D0.logpdf(W_d);
            L1_f0 = D0.logpdf(W_f);
            
            L1_s1 = D1.logpdf(W_s);
            L1_d1 = D1.logpdf(W_d);
            L1_f1 = D1.logpdf(W_f);
            
            assert(isequal(size(L1_s0), [1, n]));
            assert(isequal(size(L1_d0), [1, n]));
            assert(isequal(size(L1_f0), [1, n]));
            
            assert(isequal(L1_s0, L1_s1));
            assert(isequal(L1_d0, L1_d1));
            assert(isequal(L1_f0, L1_f1));
            
            q = deg / 2;
            lpc = q*d*log(2) + q * pdmat_lndet(S) + mvgammaln(d, q);            
            L0_s = tsuite_wishartd.my_logpdf(Sm, deg, W_s, lpc);
            L0_d = tsuite_wishartd.my_logpdf(Sm, deg, W_d, lpc);
            L0_f = tsuite_wishartd.my_logpdf(Sm, deg, W_f, lpc);
            
            devcheck('loglik (s)', L1_s0, L0_s, 1e-12);
            devcheck('loglik (d)', L1_d0, L0_d, 1e-12);
            devcheck('loglik (f)', L1_f0, L0_f, 1e-12);
        end
        
        
        function do_test_sampling(d, ty, deg)
            
            S = rand_pdmat(ty, d, 1, [0.5, 1.5]);      
            Sm = pdmat_fullform(S);
            D0 = wishartd(S, deg);
            
            N0 = 50;
            N = 2e5;
            
            Y = wishart_sample(d, deg);
            assert(isequal(size(Y), [d d]));
            assert(isequal(Y, Y'));
            assert(all(eig(Y)) > 0);
            
            W = D0.sample();
            assert(is_pdmat(W) && W.d == d && W.n == 1);
            assert(W.ty == 'f');
            assert(isequal(W.v, W.v'));
            assert(all(eig(W.v)) > 0);
            
            Ys = wishart_sample(d, deg, N);
            assert(isequal(size(Ys), [d d N]));
            for i = 1 : N0
                Yi = Ys(:,:,i);
                assert(isequal(Yi, Yi'));
                assert(all(eig(Yi)) > 0);
            end
            
            Ws = D0.sample(N);
            assert(is_pdmat(Ws) && Ws.d == d && Ws.n == N);
            for i = 1 : N0
                Wi = pdmat_fullform(Ws, i);
                assert(isequal(Wi, Wi'));
                assert(all(eig(Wi)) > 0);
            end
            
            Ymean0 = deg * eye(d);
            Wmean0 = deg * Sm;
            
            Id = eye(d);
            Yvar0 = zeros(d, d);
            Wvar0 = zeros(d, d);
            for i = 1 : d
                for j = 1 : d
                    Yvar0(i, j) = deg * (Id(i,j)^2 + Id(i,i) * Id(j,j));
                    Wvar0(i, j) = deg * (Sm(i,j)^2 + Sm(i,i) * Sm(j,j));
                end
            end
            
            Ymean = mean(Ys, 3);
            Yvar = reshape(vecvar(reshape(Ys, d*d, N)), d, d);
            Wmean = mean(Ws.v, 3);
            Wvar = reshape(vecvar(reshape(Ws.v, d*d, N)), d, d);
            
            devcheck('sample mean (std)', Ymean, Ymean0, 5e-2);
            devcheck('sample var (std)', Yvar, Yvar0, 0.3);
            devcheck('sample mean', Wmean, Wmean0, 5e-2);
            devcheck('sample var', Wvar, Wvar0, 0.3);            
        end
        
    end
    
    
    %% Auxiliary functions
    
    methods(Access='private')
        
        function run_multi(obj, tfunc)
            % run multiple test under different settings
            
            for it = numel(obj.types)
                ty = obj.types{it};
                
                for d = obj.dims
                    tfunc(d, ty, obj.deg);
                end
            end            
        end
        
    end
        
    
    methods(Static, Access='private')
                                
        function L = my_logpdf(Sm, m, W, lpc)
            
            d = W.d;
            n = W.n;
            L = zeros(1, n);
            
            for i = 1 : n
                Wi = pdmat_fullform(W, i);
                
                u1 = (m - d - 1) / 2 * lndet(Wi);
                u2 = (-1/2) * trace(Sm \ Wi);
                L(i) = u1 + u2 - lpc;
            end
        end    
    end
    
end
    
