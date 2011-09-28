classdef gaussgm < genmodel_base
    % The class that implements a Gaussian generative model
    %
    %   The formulation of the model is as follows:
    %
    %       u ~ Gauss(mu0, Cu);
    %       x ~ Gauss(A * u, Cx)
    %
    %   Here, A can be
    %       - empty:    identity transform, i.e. x ~ Gauss(u, Cx);
    %       - a scalar:
    %       - a xdim x udim matrix.
    %
    
    % Created by Dahua Lin, on Sep 27, 2011
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
            %       dx dimensional space. 
            %
            %       Here,
            %       - dx:       the dimension of observation space
            %       - Cx:       the conditional covariance of x
            %       - A:        the transformation matrix or scalar.
            %                   (If empty or omitted, using identity
            %                    transform).
            %
            
            % verify inputs
            
            if ~(isnumeric(dx) && isscalar(dx) && dx == fix(dx) && dx >= 1)
                error('gaussgm:invalidarg', ...
                    'dx should be a positive integer scalar.');
            end            
            
            if ~(is_pdmat(Cx) && Cx.d == dx && Cx.n == 1)
                error('gaussgm:invalidarg', ...
                    'Cx should be a pdmat struct with Cx.d == dx and Cx.n == 1.');
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
            model.Jx = pdmat_inv(Cx);
            model.A = A;                            
        end
        
    end
    
    
    %% Interface methods
    
    methods
           
        function tf = is_supported_optype(model, optype) %#ok<MANU>
            % Tests whether a particular operation type is supported
            %
            % This class supports 'sample', 'varinfer', and 'atom'.
            %
            
            tf = ismember(optype, {'sample', 'varinfer', 'atom'});
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
                        
            tf = isempty(pri) || (isa(pri, 'gaussd') && ...
                pri.dim == model.udim && ...
                pri.num == 1 && ...
                pri.has_ip);            
        end
        
        
        function U = posterior_params(model, pri, X, Z, optype)
            % Estimates/samples the parameters conditioned on observations
            %
            %   U = model.posterior_params(pri, X, I, 'atom');
            %       draws an atom given a subset of data points selected 
            %       by I.
            %
            %   U = model.posterior_params(pri, X, Z, 'sample');
            %       draws K parameters given grouped/weighted data set.
            %
            %   U = model.posterior_params(pri, X, Z, 'varinfer');
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
            %       - U:        a udim x 1 column vector or a udim x K
            %                   matrix.
            %
            
            Jx_ = model.Jx;
            A_ = model.A;
            
            switch optype
                case 'atom'
                    I = Z;
                    [h, J] = gaussgm_pos(pri, X(:,I), Jx_, A_, []);
                    U = gsample(h, J, 1, 'ip');
                case 'sample'                
                    [h, J] = gaussgm_pos(pri, X, Jx_, A_, Z);
                    K = size(h, 2);
                    U = zeros(size(h));
                    for k = 1 : K
                        U(:,k) = gsample(h(:,k), pdmat_pick(J, k), 1, 'ip');
                    end
                case 'varinfer'
                    [~, ~, U] = gaussgm_pos(pri, X, Jx_, A_, Z);
            end
        end
        
        
        function Lliks = evaluate_logliks(model, U, X, I)
            % Evaluates the logarithm of likelihood of samples
            %
            %   Lliks = model.evaluate_logliks(U, X)
            %
            %       evaluates the log-likelihoods of the observations
            %       in X with respect to the parameters given in U.
            %
            %       If U has K parameters, Lliks is a K x n matrix.
            %
            %   Lliks = model.evaluate_logliks(gpri, X);
            %
            %       evaluates the log marginal likelihoods with respect
            %       to the Gaussian prior given by gpri.
            %
            %   Lliks = model.evaluate_logliks(.., X, I);
            %
            %       evaluates the log likelihood values for a subset
            %       of observations selected by I.
            %
            
            if isnumeric(U)
                g = gaussd.from_mp(U, model.Cx, 'ip');                
            elseif isa(U, 'gaussd')
                g0 = U;
                g = gaussd.from_mp(g0.mu, pdmat_plus(g0.C, model.Cx), 'ip');                
            end
            
            if nargin < 4                
                Lliks = g.logpdf(X);
            else
                Lliks = g.logpdf(X(:,I));
            end
        end
        
                
        function Lpri = evaluate_logpri(model, pri, U) %#ok<MANU>
            % Evaluate the log-prior of a given set of parameters
            %
            %   Lpri = model.evaluate_logpri(pri, U);
            %
            %       evaluates the log-prior of the parameters
            %       given by U.
            %
            %       If U has K parameters, Lpri is a 1 x K vector.
            %
            
            if isempty(pri)
                Lpri = zeros(1, size(U, 2));
            else
                Lpri = pri.logpdf(U);
            end        
        end
    end
    
    
end
