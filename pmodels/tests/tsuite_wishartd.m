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
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_wishartd.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_wishartd.do_test_sampling);
        end
    end
    
    
    %% Test implementation
    
    methods(Static, Access='private')
                
        function do_test_evaluation(d, ty, deg)
            
            S = rand_pdmat(ty, d, 1, [0.5, 1.5]);
            Sm = pdmat_fullform(S);
            J = pdmat_inv(S);            
            
            n = 100;
            W_s = rand_pdmat('s', d, n, [0.5 1.5]);
            W_d = rand_pdmat('d', d, n, [0.5 1.5]);
            W_f = rand_pdmat('f', d, n, [0.5 1.5]);
            
            Ls0 = tsuite_wishartd.calc_logpdf(Sm, deg, W_s);
            Ld0 = tsuite_wishartd.calc_logpdf(Sm, deg, W_d);
            Lf0 = tsuite_wishartd.calc_logpdf(Sm, deg, W_f);
            
            Ls1 = wishartd_logpdf(S, deg, W_s);
            Ld1 = wishartd_logpdf(S, deg, W_d);
            Lf1 = wishartd_logpdf(S, deg, W_f);
            
            Ls2 = wishartd_logpdf(J, deg, W_s, [], 'inv');
            Ld2 = wishartd_logpdf(J, deg, W_d, [], 'inv');
            Lf2 = wishartd_logpdf(J, deg, W_f, [], 'inv');
            
            assert(isequal(size(Ls0), [1, n]));
            assert(isequal(size(Ld0), [1, n]));
            assert(isequal(size(Lf0), [1, n]));
                                    
            assert(isequal(size(Ls1), [1, n]));
            assert(isequal(size(Ld1), [1, n]));
            assert(isequal(size(Lf1), [1, n]));           
            
            assert(isequal(Ls1, Ls2));
            assert(isequal(Ld1, Ld2));
            assert(isequal(Lf1, Lf2));            
            
            devcheck('loglik (s)', Ls0, Ls1, 1e-12);
            devcheck('loglik (d)', Ld0, Ld1, 1e-12);
            devcheck('loglik (f)', Lf0, Lf1, 1e-12);
        end
        
        
        function do_test_sampling(d, ty, deg)
            
            S = rand_pdmat(ty, d, 1, [0.5, 1.5]);      
            Sm = pdmat_fullform(S);
            
            N0 = 50;
            N = 2e5;           
            
            Ys = wishartd_sample(d, deg, N);
            assert(isequal(size(Ys), [d d N]));
            for i = 1 : N0
                Yi = Ys(:,:,i);
                assert(isequal(Yi, Yi'));
                assert(all(eig(Yi)) > 0);
            end
            
            Ws = wishartd_sample(S, deg, N);
            assert(isequal(size(Ws), [d d N]));
            for i = 1 : N0
                Wi = Ws(:,:,i);
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
            Wmean = mean(Ws, 3);
            Wvar = reshape(vecvar(reshape(Ws, d*d, N)), d, d);
            
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
                                
        function L = calc_logpdf(Sm, df, W)
            
            d = W.d;
            lpc = df*d*log(2)/2 + df * lndet(Sm)/2 + mvgammaln(d, df/2);             
            n = W.n;
            L = zeros(1, n);
            
            for i = 1 : n
                Wi = pdmat_fullform(W, i);
                
                u1 = (df - d - 1) / 2 * lndet(Wi);
                u2 = (-1/2) * trace(Sm \ Wi);
                L(i) = u1 + u2 - lpc;
            end
        end    
    end
    
end
    
