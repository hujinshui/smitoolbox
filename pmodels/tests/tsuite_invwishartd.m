classdef tsuite_invwishartd
    % Test suite for invwishartd
    %
    
    % Created by Dahua Lin, on Sep 3, 2011

    %% Properties
    
    properties
        dims = [1 2 3 5];
        types = {'s', 'd', 'f'};
        deg = 8.5;        
    end
        
    %% Test cases
    
    methods
    
        function test_basics(obj)
            run_multi(obj, @tsuite_invwishartd.do_test_basics);
        end
        
        function test_statistics(obj)
            run_multi(obj, @tsuite_invwishartd.do_test_statistics);
        end
        
        function test_evaluation(obj)
            run_multi(obj, @tsuite_invwishartd.do_test_evaluation);
        end
        
        function test_sampling(obj)
            run_multi(obj, @tsuite_invwishartd.do_test_sampling);
        end
    end
    
    
    %% Test implementation
    
    methods(Static, Access='private')
    
        function do_test_basics(d, ty, deg)
                        
            Phi = rand_pdmat(ty, d, 1, [0.5, 1.5]);
            
            D0 = invwishartd(Phi, deg);
            D1 = invwishartd(Phi, deg, 'pre');
            
            assert(D0.num == 1);
            assert(D0.dim == d);
            assert(D0.deg == deg);
            assert(isequal(D0.Phi, Phi));
            assert(isempty(D0.invPhi));
            assert(isempty(D0.lpconst));
            
            assert(D1.num == 1);
            assert(D1.dim == d);
            assert(D1.deg == deg);
            assert(isequal(D1.Phi, Phi));
            assert(~isempty(D1.invPhi));
            assert(~isempty(D1.lpconst));
            
            % verify pre-computation
            
            q = deg / 2;
            lpc = q*d*log(2) - q * pdmat_lndet(Phi) + mvgammaln(d, q);
            
            devcheck('lpconst', D1.lpconst, lpc, 1e-13);
            assert(isequal(pdmat_inv(Phi), D1.invPhi));
        end
        
        
        function do_test_statistics(d, ty, deg)
                                    
            Phi = rand_pdmat(ty, d, 1, [0.5, 1.5]);
            Phim = pdmat_fullform(Phi);
            
            D0 = invwishartd(Phi, deg);
            
            Mean0 = Phim / (deg - d - 1);
            Mode0 = Phim / (deg + d + 1);
            
            devcheck('mean (D0)', pdmat_fullform(mean(D0)), Mean0, 1e-13);
            devcheck('mode (D0)', pdmat_fullform(mode(D0)), Mode0, 1e-13);
        end
    
        
        function do_test_evaluation(d, ty, deg)
            
            Phi = rand_pdmat(ty, d, 1, [0.5, 1.5]);
            Phim = pdmat_fullform(Phi);
            
            D0 = invwishartd(Phi, deg);
            D1 = invwishartd(Phi, deg, 'pre');            
            
            n = 100;
            C_s = rand_pdmat('s', d, n, [0.5, 1.5]);
            C_d = rand_pdmat('d', d, n, [0.5, 1.5]);
            C_f = rand_pdmat('f', d, n, [0.5, 1.5]);
            
            L1_s0 = D0.logpdf(C_s);
            L1_d0 = D0.logpdf(C_d);
            L1_f0 = D0.logpdf(C_f);
            
            L1_s1 = D1.logpdf(C_s);
            L1_d1 = D1.logpdf(C_d);
            L1_f1 = D1.logpdf(C_f);
            
            assert(isequal(size(L1_s0), [1, n]));
            assert(isequal(size(L1_d0), [1, n]));
            assert(isequal(size(L1_f0), [1, n]));
            
            assert(isequal(L1_s0, L1_s1));
            assert(isequal(L1_d0, L1_d1));
            assert(isequal(L1_f0, L1_f1));
                                    
            q = deg / 2;
            lpc = q*d*log(2) - q * pdmat_lndet(Phi) + mvgammaln(d, q);
            L0_s = tsuite_invwishartd.my_logpdf(Phim, deg, C_s, lpc);
            L0_d = tsuite_invwishartd.my_logpdf(Phim, deg, C_d, lpc);
            L0_f = tsuite_invwishartd.my_logpdf(Phim, deg, C_f, lpc);
            
            devcheck('loglik (s)', L1_s0, L0_s, 1e-12);
            devcheck('loglik (d)', L1_d0, L0_d, 1e-12);
            devcheck('loglik (f)', L1_f0, L0_f, 1e-12);
        end
        
        
        function do_test_sampling(d, ty, deg)
                                    
            Phi = rand_pdmat(ty, d, 1, [0.5, 1.5]);            
            D0 = invwishartd(Phi, deg);
            
            N0 = 50;
            N = 2e5;
            
            C = D0.sample();
            assert(is_pdmat(C) && C.d == d && C.n == 1);
            assert(C.ty == 'f');
            assert(isequal(C.v, C.v'));
            assert(all(eig(C.v)) > 0);
            
            Cs = D0.sample(N);
            assert(is_pdmat(Cs) && Cs.d == d && Cs.n == N);
            for i = 1 : N0
                Ci = pdmat_fullform(Cs, i);
                assert(isequal(Ci, Ci'));
                assert(all(eig(Ci)) > 0);
            end
            
            Cmean0 = pdmat_fullform(mean(D0));
            Cmean = mean(Cs.v, 3);
            devcheck('sample mean', Cmean, Cmean0, 5e-3);
            
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
        
        function L = my_logpdf(Phim, m, C, lpc)
            
            d = C.d;
            n = C.n;
            L = zeros(1, n);
            
            for i = 1 : n
                Ci = pdmat_fullform(C, i);
                
                u1 = - (m + d + 1) / 2 * lndet(Ci);
                u2 = (-1/2) * trace(Ci \ Phim);
                L(i) = u1 + u2 - lpc;
            end
        end
    end
    
    
end

