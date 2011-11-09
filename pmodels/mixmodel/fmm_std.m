classdef fmm_std < smi_prg
    % The class that implements standard finite mixture model     
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Sep 27, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')
        gm;         % The underlying generative model
        pri;        % The parameter prior 
        K;          % The number of mixture components
        jak = 0;    % whether to use jacket GPU computation
    end
    
    
    properties
        dalpha = 1; % the concentration parameter of the Dirichlet prior of pi
                    % dird_alpha >= 1
                    % one can set dird_alpha to inf, in that case
                    % pi is fixed to [1/k, ..., 1/k].        
    end
    
    methods
        function prg = set.dalpha(prg, v)
            if ~(isfloat(v) && isreal(v) && v >= 1)
                error('fmm_std:invalidarg', ...
                    'The value of dalpha should be a real numberin [1, inf]');
            end
            prg.dalpha = v;
        end
    end
    
    
    %% Construction
    
    methods
        
        function prg = fmm_std(gm, pri, K, op)
            % Creates a standard Gaussian mixture model (GMM) program
            %                        
            %   prg = gmm_std(gm, pri, K);
            %
            %       Construsts an SMI program that implements the 
            %       inference of FMM parameters.
            %
            %       Input arguments:
            %       - gm:       the underlying generative model, which
            %                   should be an instance of a derived class
            %                   of genmodel_base.
            %
            %       - pri:      the prior distribution of the parameters,
            %                   which should be in a form acceptable by
            %                   gm.
            %
            %       - K:        the number of mixture components to be
            %                   estimated.
            %
            %   prg = gmm_std(gm, pri, K, 'jak');
            %
            %       Tells the program to use Jacket GPU computation.
            %       
                        
            % verify input arguments
            
            if ~isa(gm, 'genmodel_base')
                error('fmm_std:invalidarg', ...
                    'gm should be an instance of a derived class of genmodel_base.');
            end
            
            if ~gm.is_valid_prior(pri)
                error('fmm_std:invalidarg', 'The input pri is not valid.');
            end
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
                error('fmm_std:invalidarg', ...
                    'K should be a positive integer scalar.');
            end
            
            if nargin >= 4 && strcmpi(op, 'jak')
                prg.jak = 1;
            end
            
            % create object
            
            prg.gm = gm;
            prg.pri = pri;
            prg.K = double(K);
        end
        
    end
    
    
    %% Implementation
    
    methods
    
        function [Sd, Sc] = initialize(prg, obs, L0, optype)
            % Initialize the states
            %
            %   [Sd, Sc] = prg.initialize(obs, L0, 'sample');
            %   [Sd, Sc] = prg.initialize(obs, L0, 'varinfer');
            %
            %       Input arguments:
            %       - obs:  The observation set (in a form acceptable by
            %               gm).
            %
            %       - L0:   It can be given in three forms:
            %               - empty:  use default method to initialize
            %               - 1 x n label vector
            %               - K x n label distribution matrix            
            %
            %       Note that these two arguments are packed in a 
            %       1 x 2 cell array as a single input.
            %
            %       Output arguments:
            %       - Sd:   dynamic state struct, with fields: 
            %               - params:   the parameters
            %               - aux:      the auxiliary structure created
            %                           by the underlying model
            %               - Liks:     K x n pairwise likelihood table
            %               - Pi:       prior probabilities of components
            %                           (1 x K vector, initialized to
            %                            a uniform distribution)
            %
            %           If optype is 'sample', Sd also has
            %               - Z:    1 x n label vector
            %               - grps: 1 x K cell array of grouped indices
            %
            %           If optype is 'varinfer', Sd also has
            %               - Q:    K x n matrix: posterior label
            %                       distributions.
            %               
            %       - Sc:   fixed state struct, with fields:
            %               - obs:      the observation set
            %               - n:        the number of samples in obs
            %               - K:        the number of mixture components
            %               - optype;   the operation type string
            %               - do_samp;  whether it does sampling            
            %
            
            % verify inputs
            
            gmdl = prg.gm;
            K_ = prg.K;
            
            n = gmdl.get_num_observations(obs);                                  
            
            if ~isempty(L0)
                if ~(isnumeric(L0) && ndims(L0) == 2 && ...
                        size(L0,2) == n && (size(L0,1) == 1 || size(L0,1) == K_));
                    error('fmm_std:invalidarg', ...
                        'L0 should be a 1 x n row vector or K x n matrix.');                
                end
                
                if size(L0, 1) == 1
                    if ~all(L0 == fix(L0) && L0 >= 1 && L0 <= K_)
                        error('fmm_std:invalidarg', ...
                            'Some labels in L0 is invalid.');
                    end
                end                                        
            end
                            
            switch optype
                case 'sample'                    
                    do_samp = true;
                case 'varinfer'
                    do_samp = false;
                otherwise
                    error('fmm_std:invalidarg', ...
                        'Unsupported operation type %s', optype);
            end
            
            if ~gmdl.is_supported_optype(optype)
                error('fmm_std:invalidarg', ...
                    'The optype %s is not supported by the underlying model', optype);
            end
            
            if prg.jak
                if ~strcmp(optype, 'varinfer')
                    error('fmm_std:invalidarg', ...
                        'Only varinfer can be used with jacket GPU.');
                end
            end
                                
            % create states
            
            Sc.obs = obs;
            Sc.n = n;
            Sc.K = prg.K;
            Sc.optype = optype;
            Sc.do_samp = do_samp; 
            
            Sd.params = []; 
            Sd.aux = [];
            Sd.Liks = [];
            Sd.Pi = constmat(1, Sc.K, 1 / Sc.K);
                        
            if do_samp
                if isempty(L0)
                    Sd.Z = rand_label(Sc.K, n, 'v');
                elseif size(L0, 1) == 1
                    Sd.Z = L0;
                else
                    [~, Sd.Z] = max(L0, [], 1);
                end
                Sd.grps = intgroup(Sc.K, Sd.Z);
            else
                if isempty(L0)
                    Sd.Q = rand_label(Sc.K, n, 'b');
                elseif size(L0, 1) == 1
                    Sd.Q = l2mat(Sc.K, L0);
                else
                    Sd.Q = L0;
                end
            end 
            
            if prg.jak
                Sd.Pi = gdouble(Sd.Pi);
                Sd.Q = gdouble(Sd.Q);                 
            end
            
        end    
        
        
        function [Sd, Sc] = update(prg, Sd, Sc)
            % Updates the states (one iteration)
            %
            %   [Sd, Sc] = prg.update(Sd, Sc);
            %
                        
            gmdl = prg.gm;
            prior = prg.pri;
            X = Sc.obs; 
            K_ = Sc.K;
            samp = Sc.do_samp;

            % estimate component parameters
            
            if samp
                Z0 = Sd.grps;                
            else
                Z0 = Sd.Q;
            end
            
            [Sd.params, Sd.aux] = gmdl.posterior_params(prior, X, Sd.aux, Z0, Sc.optype);
            Sd.Liks = gmdl.evaluate_logliks(Sd.params, X, Sd.aux);
            
            % infer labels                                    
            
            Q = fmm_inferQ(Sd.Pi.', Sd.Liks);
            
            if samp
                Sd.Z = ddsample(Q, 1);
                Sd.grps = intgroup(K_, Sd.Z);
            else
                Sd.Q = Q;
            end
            
            % estimate Pi
            
            if samp
                Sd.Pi = fmm_estPi(K_, Sd.Z, prg.dalpha, 'sample').';
            else
                Sd.Pi = fmm_estPi(K_, Q, prg.dalpha).';
            end            
            
        end
          
                
        function Sp = make_output(prg, Sd, Sc) %#ok<MANU>
            % Makes output struct from given states
            %
            %   Sp = prg.make_output(Sd, Sc);
            %       
            %       The output struct has the following fields:
            %       - params:   the estimated/sampled component params
            %       - Pi:       the estimated/sampled prior probability of
            %                   components. [1 x K]
            %
            %       If it is doing sampling, Sp also has:
            %       - Z:    the label vector [1 x n]
            %
            %       If it is doing variational EM, Sp also has
            %       - Q:    the posterior distribution of labels [K x n]
            %
            
            Sp.params = Sd.params;            
            Sp.Pi = Sd.Pi;
            
            if Sc.do_samp
                Sp.Z = Sd.Z;
            else
                Sp.Q = Sd.Q;
            end                               
        end
        
        
        function objv = evaluate_objective(prg, Sd, Sc)
            % Evaluates the objective value w.r.t the given states
            %
            %   objv = prg.evaluate_objective(Sd, Sc);
            %
            
            gmdl = prg.gm;
            prior = prg.pri;
            K_ = Sc.K;
            n = Sc.n;
            
            % log-prior of component params
            
            lpri_u = gmdl.evaluate_logpri(prior, Sd.params, Sd.aux);
            lpri_u = double(lpri_u);
            
            % log-prior of labeling and Pi
            
            logPi = log(Sd.Pi);
            
            if Sc.do_samp
                tw = intcount(K_, Sd.Z).';
            else
                tw = sum(Sd.Q, 2);
            end
            lpri_z = logPi * tw;  
            lpri_z = double(lpri_z);
            
            da = prg.dalpha;
            if isfinite(da) && da > 1
                lpri_pi = sum(logPi) * (da - 1);
                lpri_pi = double(lpri_pi);
            else
                lpri_pi = 0;
            end
            
            % log-likelihood
            
            if Sc.do_samp                
                llik = Sd.Liks(sub2ind([K_, n], Sd.Z, 1:n));
                llik = sum(llik);                
            else
                llik = dot(Sd.Liks, Sd.Q, 1);
                llik = sum(llik);
            end      
            llik = double(llik);
            
            % labeling entropy
            
            if Sc.do_samp
                ent_z = 0;
            else
                ent_z = sum(ddentropy(Sd.Q));
                ent_z = double(ent_z);
            end
            
            % combine
            
            objv = lpri_u + lpri_z + lpri_pi + llik + ent_z;
            
        end
            
    end
    
end




