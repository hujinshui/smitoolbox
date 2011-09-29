classdef gaussgm < genmodel_base
    % The class that implements a Gaussian generative model
    %
    %   The formulation of the model is as follows:
    %
    %     Formulation 1
    %
    %       u ~ Gauss(mu0, Cu);  
    %       x ~ Gauss(A * u, Cx); 
    %
    %     Formulation 2
    %
    %       u ~ Gauss(mu0, Cu);
    %       Cx ~ ScaledInverseChiSquare or InverseWishart
    %       x ~ Gauss(A * u, Cx);
    %
    %   Here, A can be
    %       - empty:    identity transform, i.e. x ~ Gauss(u, Cx);
    %       - a scalar:
    %       - a xdim x udim matrix.
    %
    %
    %   The prior distribution should be given in the following form:
    %   empty or gpri for Formulation 1, or {gpri, Cpri} for Formulation 2.
    %
    %   - gpri:     The prior distribution of the mean vector u.
    %               It is either empty (uninformative prior), or 
    %               a Gaussian distribution with gpri.dim == udim,
    %               gpri.num == 1, and gpri.has_ip.
    %
    %   - Cpri:     The prior distribution of the covariance matrix.
    %   
    %               - a char specifying the form of Cx to be estimated, 
    %                 which can be either 's', 'd', or 'f'. This indicates
    %                 an uninformative prior of the corresponding form.
    %
    %               - an object of class scale_invchi2d, as an scaled 
    %                 inverse chi-squared prior. Specifically,
    %                   if Cpri.d == 1, Cx will be in scale form, or 
    %                   if Cpri.d == udim, Cx will be in diagonal form.
    %                   
    %               - an object of class invwishartd, as an inversed 
    %                 Wishart distribution prior. In this case, Cx will
    %                 be in full matrix form.
    %   
    %
    
    % History
    % -------
    %   - Created by Dahua Lin, on Sep 27, 2011
    %   - Modified by Dahua Lin, on Sep 28, 2011
    %       - add the support of Cx estimation
    %
    
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')
        
        udim;   % the dimension of u
        xdim;   % the dimension of x
        
        Cx;     % the covariance of the conditional Gaussian distribution
        Jx;     % the inverse of Cx
        
        A;      % the transformation matrix
                % A can be empty, indicating identity transform
                % A can also be a scalar, indicating scaling transform
                
        tied_cov;   % whether the covariance of different components
                    % are tied to the same
    end        
    
    
    
    %% Construction
    
    methods
        
        function model = gaussgm(dx, Cx, A)
            % Constructs a Gaussian generative model
            %                    
            %   model = gaussgm(dx, Cx);
            %   model = gaussgm(dx, Cx, A);
            %
            %       constructs a Gaussian generative model over a 
            %       dx dimensional space, based on Formulation 1.
            %
            %   model = gaussgm(dx);
            %   model = gaussgm(dx, [], A);
            %
            %       constructs a Gaussian generative model over a 
            %       dx dimensional space, based on Formulation 2.
            %
            %   model = gaussgm(dx, 'tied_cov');
            %   model = gaussgm(dx, 'tied_cov', A);
            %
            %       constructs a Gaussian generative model over a 
            %       dx dimensional space, based on Formulation 2,
            %       and with the setting that the covariance is
            %       shared by different components.            
            %
            %       Here,
            %       - dx:       the dimension of observation space
            %       - Cx:       the conditional covariance of x
            %       - A:        the transformation matrix or scalar.
            %                   (If empty or omitted, using identity
            %                    transform).
            %
            %       When Cx is given, the model uses Formulation 1,
            %       otherwise, when Cx is omitted or input as empty,
            %       the model uses Formulation 2.
            %
            
            % verify inputs
            
            if ~(isnumeric(dx) && isscalar(dx) && dx == fix(dx) && dx >= 1)
                error('gaussgm:invalidarg', ...
                    'dx should be a positive integer scalar.');
            end            
            
            if nargin < 2 || isempty(Cx)
                Cx = [];
                c_tied = false;
                
            elseif ischar(Cx)
                if ~strcmpi(Cx, 'tied_cov')
                    error('gaussgm:invalidarg', ...
                        'The 2nd argument in constructing gaussgm is invalid.');
                end
                Cx = [];
                c_tied = true;
                
            else
                if ~(is_pdmat(Cx) && Cx.d == dx && Cx.n == 1)
                    error('gaussgm:invalidarg', ...
                        'Cx should be a pdmat struct with Cx.d == dx and Cx.n == 1.');
                end
                c_tied = true;
            end
                        
            if nargin >= 3 && ~isempty(A)
                if ~(isfloat(A) && isreal(A) && (isscalar(A) || ...
                        (ndims(A) == 2 && size(A,1) == dx)) )
                    error('gaussgm:invalidarg', ...
                        'A should be either a scalar or a real matrix with dx rows.');
                end
                if isscalar(A)
                    du = dx;
                else
                    du = size(A, 2);
                end
            else
                A = [];
                du = dx;
            end
               
            % set fields
            
            model.udim = du;
            model.xdim = dx;
            model.Cx = Cx;
            if ~isempty(Cx)
                model.Jx = pdmat_inv(Cx);
            end
            model.A = A;   
            model.tied_cov = c_tied;
        end
        
    end
    
    
    %% Interface methods
    
    methods
           
        function tf = is_supported_optype(model, optype) 
            % Tests whether a particular operation type is supported
            %
            % This class supports 'sample', 'varinfer', and 'atom'.
            %
            
            if ~isempty(model.Cx)
                tf = ismember(optype, {'sample', 'varinfer', 'atom'});
            else
                tf = ismember(optype, {'sample', 'varinfer'});
            end
        end
        
        
        function n = get_num_observations(model, X)
            % Gets the number of samples in the observation set
            %
            %   n = model.get_num_observations(X);
            %       verifies the validity of X as an observation set,
            %       and returns the number of samples in X,
            %
            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == model.xdim)
                error('gaussgm:invalidobs', ...
                    'X is not a valid observation set.');
            end            
            n = size(X, 2);
        end
        
        
        function tf = is_valid_prior(model, pri)
            % Tests the validity of a prior
            %
            %   tf = is_valid_prior(model, pri);
            %
            %       it returns true when pri is a valid prior, i.e.
            %       pri is an instance of gaussd, with pri.dim == udim,
            %       and pri.num == 1. In addition, pri.has_ip must be
            %       true.
            %
            %       pri = [] is also valid, which indicates an
            %       uninformative prior.
            %
                        
            tf = gaussgm.verify_prior(model.udim, isempty(model.Cx), pri);
        end
        
        
        function [params, aux] = posterior_params(model, pri, X, aux, Z, optype)
            % Estimates/samples the parameters conditioned on observations
            %
            %   [params, aux] = model.posterior_params(pri, X, aux, I, 'atom');
            %       draws an atom given a subset of data points selected 
            %       by I.
            %
            %   [params, aux] = model.posterior_params(pri, X, aux, Z, 'sample');
            %       draws K parameters given grouped/weighted data set.
            %
            %   [params, aux] = model.posterior_params(pri, X, aux, Z, 'varinfer');
            %       estimates/infer K parameters given grouped/weighted 
            %       data set.
            %
            %       Input Arguments:
            %       - pri:      the Gaussian prior of the parameters
            %       - X:        the observation data set
            %       - I:        an index vector.
            %       - Z:        It can be either a cell array of grouped
            %                   indices, or a K x n matrix of weights.
            %
            %       Output Arguments:
            %
            %       For Formulation 1
            %       - params:   a udim x K matrix.
            %
            %       For Formulation 2
            %       - params:   a struct with the following fields:
            %                   'U':    the mean vector(s) [udim x K]
            %                   'Cx':   the estimated covariance,
            %                           in form of pdmat struct.
            %
            %       For both formulation:
            %
            %       - aux:      a struct with the following fields:
            %                   'upri':     the prior of u
            %                   'cpri':     the prior of Cx
            %                   'cf':       the form of Cx to be estimated
            %                   'Cx':       currently estimated Cx
            %                               (only for Formulation 2)
            %
            
            % basic info
            
            du = model.udim;
            estCx = isempty(model.Cx);                            
            
            % prepare aux
            
            if isempty(aux) % first time invocation
                [tf, upri, cpri, cf] = gaussgm.verify_prior(du, estCx, pri);
                if ~tf
                    error('gaussgm:rterror', 'pri is invalid.');
                end
                aux = [];
                aux.upri = upri;
                aux.cpri = cpri;
                aux.cf = cf;
                aux.Cx = [];
            else
                upri = aux.upri;
                cpri = aux.cpri;
                cf = aux.cf;
                Cx_ = aux.Cx;
            end
            
            % Get Jx_
            
            if estCx;
                if isempty(Cx_)
                    Cx_ = gaussgm.bootCx(cf, cpri);
                end
                Jx_ = pdmat_inv(Cx_);
            else
                Jx_ = model.Jx;
            end
            
            % Estimate mean vectors
            
            A_ = model.A;                        
            
            switch optype
                case 'atom'
                    I = Z;
                    [h, J] = gaussgm_pos(upri, X(:,I), Jx_, A_, []);
                    U = gsample(h, J, 1, 'ip');
                case 'sample'                
                    [h, J] = gaussgm_pos(upri, X, Jx_, A_, Z);
                    K = size(h, 2);
                    U = zeros(size(h));
                    for k = 1 : K
                        U(:,k) = gsample(h(:,k), pdmat_pick(J, k), 1, 'ip');
                    end
                case 'varinfer'
                    [~, ~, U] = gaussgm_pos(upri, X, Jx_, A_, Z);
            end
            
            % Estimate covariance matrix
            
            if estCx
                switch optype
                    case 'atom'
                        [~, ~, Cx_e] = ...
                            gaussgm_cov_pos(cf, cpri, X(:,I), U, [], c_tied, 'sample');
                    case 'sample'
                        [~, ~, Cx_e] = ...
                            gaussgm_cov_pos(cf, cpri, X, U, Z, c_tied, 'sample');
                    case 'varinfer'
                        [~, ~, Cx_e] = ...
                            gaussgm_cov_pos(cf, cpri, X, U, Z, c_tied);
                end                
            end
            
            % Form output
            
            if ~estCx
                params = U;
            else
                params = [];
                params.U = U;
                params.Cx = Cx_e;
                aux.Cx = Cx_e;
            end
            
        end
        
        
        function Lliks = evaluate_logliks(model, params, X, ~, I)
            % Evaluates the logarithm of likelihood of samples
            %
            %   Lliks = model.evaluate_logliks(params, X, aux)
            %
            %       evaluates the log-likelihoods of the observations
            %       in X with respect to the parameters given in U.
            %
            %       If U has K parameters, Lliks is a K x n matrix.
            %
            %   Lliks = model.evaluate_logliks(gpri, X, aux);
            %
            %       evaluates the log marginal likelihoods with respect
            %       to the Gaussian prior given by gpri.
            %
            %   Lliks = model.evaluate_logliks(.., X, aux, I);
            %
            %       evaluates the log likelihood values for a subset
            %       of observations selected by I.
            %
            
            if isnumeric(params)
                U = params;
                g = gaussd.from_mp(U, model.Cx, 'ip');  
                
            elseif isa(params, 'gaussd')
                g0 = params;
                g = gaussd.from_mp(g0.mu, pdmat_plus(g0.C, model.Cx), 'ip');                
            end
            
            if nargin < 5              
                Lliks = g.logpdf(X);
            else
                Lliks = g.logpdf(X(:,I));
            end
        end
        
                
        function lpri = evaluate_logpri(model, pri, params, ~) %#ok<MANU>
            % Evaluate the total log-prior of a given set of parameters
            %
            %   Lpri = model.evaluate_logpri(pri, params, aux);
            %
            %       evaluates the log-prior of the parameters
            %       given by U.
            %
            %       If U has K parameters, Lpri is a 1 x K vector.
            %
            
            if isempty(pri)
                lpri = 0;
                
            else
                if isnumeric(params)
                    U = params;
                    lpri = sum(pri.logpdf(U));
                end        
            end
        end
    end
    
    
    %% Private methods
    
    methods(Static, Access='private')
        
        function [tf, upri, cpri, cf] = verify_prior(du, estCx, pri)
            % Verifies the validity of the prior
            %
            %   tf:         whether it is valid
            %   upri:       the prior object for u, or empty
            %   cpri:       the prior object for Cx, or empty
            %   cf:         the form of Cx
            %            
            
            if ~estCx   % Formualtion 1: Cx fixed
                
                if isempty(pri)
                    tf = true;
                    upri = [];
                elseif isa(pri, 'gaussd') 
                    tf = pri.dim == du && pri.num == 1 && pri.has_ip;
                    upri = pri;
                else
                    tf = false;
                    upri = [];
                end
                cpri = [];
                cf = char(0);                
                
            else    % Formulation 2: Cx to be estimated
                
                if isempty(pri)
                    tf = true;
                    upri = [];
                    cpri = [];
                    cf = 'f';
                    
                elseif iscell(pri) && numel(pri) 
                    e1 = pri{1};
                    e2 = pri{2};
                    
                    if isempty(e1)
                        tf1 = true;
                        upri = [];
                    elseif isa(e1, 'gaussd')
                        tf1 = pri.dim == du && pri.num == 1 && pri.has_ip;
                        upri = e1;
                    else
                        tf1 = false;
                        upri = [];
                    end
                    
                    if isempty(e2)
                        tf2 = true;
                        cpri = [];
                        cf = 'f';
                    elseif ischar(e2) && isscalar(e2)
                        tf2 = (e2 == 's' || e2 == 'd' || e2 == 'f');
                        cpri = [];
                        cf = e2;
                    elseif isa(e2, 'invgammad')
                        tf2 = (e2.num == 1 && (e2.dim == 1 || e2.dim == du));
                        cpri = e2;
                        if e2.dim == 1
                            cf = 's';
                        else
                            cf = 'd';
                        end
                    elseif isa(e2, 'invwishartd')
                        tf2 = e2.num == 1 && e2.dim == du;
                        cpri = e2;
                        cf = 'f';
                    else
                        tf2 = false;
                        cpri = [];
                        cf = char(0);
                    end
                    
                    tf = tf1 && tf2;
                    
                else
                    tf = false;
                    upri = [];
                    cpri = [];
                    cf = char(0);                    
                end                                    
            end
        end
    
        
        function Cx = bootCx(d, cf, cpri)
            % Bootstrap Cx (for Formulation 2)
            %
            
            switch cf
                case 's'
                    if isempty(cpri)
                        v = 1;
                    else
                        v = cpri.mode();
                    end
                    Cx = pdmat('s', d, v);
                    
                case 'd'
                    if isempty(cpri)
                        v = ones(d, 1);
                    else
                        v = cpri.mode();
                    end
                    Cx = pdmat('d', d, v);
                    
                case 'f'
                    if isempty(cpri)
                        C = eye(d);
                    else
                        C = cpri.mode();
                    end
                    Cx = pdmat('d', d, C);                    
            end            
        end
        
        
    end
             
end


