classdef kmeans_alg
    % The class to implement the K-means algorithm
    %
    % This class implementes Standard K-means algorithm in two ways:
    %
    %   - Conventional way: it computes the distances between each pair
    %     of sample and center at each iteration.
    %
    %   - Accelerated way: it avoids redundant distance computation 
    %     using triangle inequality, as described by the following paper:
    %
    %       Charles Elkan. "Using the Triangle Inequality to Accelerate 
    %       k-Means". Proceedings of 20th International Conference on
    %       Machine Learning (ICML-03), Washington DC. 2003.
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Sep 27, 2010
    %
    
    
    properties
        dist_func = @km_euc;
        
        max_iter = 300;     % the maximum number of iterations
        tol_c = 0;          % the maximum allowable changes of labels at convergence
        display = 'off';    % the level of information displaying                
    end
    
    
    methods                      
        
        function [L, M, info] = run_std(alg, X, M0)
            % Runs K-means in the conventional way
            %
            
            %% intialization
            
            cfunc = alg.cost_func;
            mfunc = alg.mean_func;
            
            M = M0;
            K = size(M, 2);
            Gp = cell(1, K);
            
            cost = cfunc(M, X);
            [mcv, L] = min(cost, [], 1);
            G = intgroup([1, K], L);
            
            
            %% iterations
            
            maxit = alg.max_iter;
            it = 0;
            converged = false;
            
            while ~converged && it < maxit 
                
                it = it + 1;
                
                % update centers                                
                
                M = kmeans_alg.update_means(M, X, Gp, G, mfunc);
                                
                % update labeling
                
                cost = cfunc(M, X);
                Lp = L;
                [mcv, L] = min(cost, [], 1);
                
                if ~isequal(Lp, L)
                    Gp = G;
                    G = intgroup([1, K], L);
                else
                    converged = true;
                end                                
            end
            
            %% output info
            
            if nargout >= 3
                info.converged = converged;
                info.niters = it;
                info.totalcost = sum(mcv);
            end

        end               
        
        
    end
    
    
    methods(Static, Access='private')        
        
        function M = update_means(M, X, Gp, G, mfunc)
            % update mean vector(s)
            %
            %   Gp: previous grouping
            %   G:  current grouping
            %
            
            K = size(M, 2);
            
            for k = 1 : K                
                gk = G{k};
                if ~isequal(Gp{k}, gk) && ~isempty(gk)                    
                    M(:, k) = mfunc(X(:, gk));   
                end                
            end
        end
        
    end
    
    
end