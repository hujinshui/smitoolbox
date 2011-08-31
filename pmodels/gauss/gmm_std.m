function prg = gmm_std(X, K, Cx, gpri, iwspri, method)
% Creates a standard Gaussian mixture model (GMM) program
%
%   The GMM program is an SMI program that implements the inference
%   over the following formulation.
%
%       u_k ~ Gauss(mu, Cu);        for k = 1, ..., K
%       x_i ~ Gauss(u_{z_i}, Cx);   for i = 1, ..., N
%
%   Cx is either given, or from an inverse Wishart distribution.
%
%   prg = gmm_std(X, K);
%   prg = gmm_std(X, K, Cx);
%   prg = gmm_std(X, K, [], gpri);
%   prg = gmm_std(X, K, Cx, gpri);
%
%   prg = gmm_std(X, K, [], gpri, iwspri); (not implemented yet)
%
%   prg = gmm_std( ......, method);
%       
%       Construsts an SMI program that implements the inference of
%       GMM parameters. 
%       
%       Input arguments:
%       - X:        the observed data [xdim x n]
%       - K:        the number of mixture components
%       - Cx:       the covariance for generating x.
%                   (If left empty or omitted, it is to be inferred)
%
%       - gpri:     the Gaussian prior for the mean
%       - iwspri:   the Inverse Wishart prior for the covariance
%
%       - method:   the inference method:
%                   - 'em': using expectation-maximization (default)
%                   -' gs': using gibbs sampling
%

% Created by Dahua Lin, on Aug 31, 2011
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2 && ~issparse(X) && isreal(X))
    error('gmm_std:invalidarg', 'X should be a real matrix.');
end
[d, N] = size(X);

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
    error('gmm_std:invalidarg', 'K should be a non-negative integer scalar.');
end

if nargin < 3 || isempty(Cx)
    est_Cx = true;
    Cx = [];
else
    est_Cx = false;
    if ~(is_pdmat(Cx) && Cx.d == d)
        error('gmm_std:invalidarg', ...
            'Cx should be a pdmat struct with Cx.dim == size(X,1).');
    end
end

if nargin < 4 || isempty(gpri)
    mu = [];
    Cu = [];
else
    if ~(isa(gpri, 'gaussd') && gpri.has_mp)
        error('gmm_std:invalidarg', ...
            'gpri should be a gaussd object with mean parameters.');
    end
    if ~(gpri.num == 1 && gpri.dim == d)
        error('gmm_std:invalidarg', 'gpri is invalid.');
    end
    mu = gpri.mu;
    Cu = gpri.C;
    
    if isequal(mu, 0)
        mu = zeros(d, 1);
    end
end
       
if nargin < 6
    method = 'em';
else
    if ~isvarname(method)
        error('gmm_std:invalidarg', 'method is not a valid name.');
    end
    method = lower(method);
    if ~(strcmp(method, 'em') || strcmp(method, 'gs'))
        error('gmm_std:invalidarg', 'Unknown method name %s', method);
    end
end

if est_Cx
    error('gmm_std:notimplement', ...
        'Estimation of Cx has not yet been implemented.');
end


%% main

% create core model

if isempty(Cx)
    if isempty(Cu)
        gm = gaussgm(d, d);
    else
        gm = gaussgm(d, Cu);
    end
else
    if isempty(Cu)
        gm = gaussgm(Cx, d);
    else
        gm = gaussgm(Cx, Cu);
    end
end

has_u_pri = ~isempty(Cu);

% create the framework

frmwork = smi_frmwork();

% add variables

frmwork.add_var('X', 'double', [d, N]);     % observations
frmwork.add_var('U', 'double', [d, K]);    % mean vectors

if has_u_pri
    frmwork.add_var('mu', 'double', [d, 1]);    
end

if est_Cx
    frmwork.add_var('Cx', 'struct', 1);     % covariance matrix
end

switch method
    case 'em'
        frmwork.add_var('Z', 'double', [K, N]);
    case 'gs'
        frmwork.add_var('Z', 'double', [1, N]);
end

% add functions

switch method
    case 'em'
        infer_method = 'mapest';
    case 'gs'
        infer_method = 'sample';
end

mean_inferrer = gaussgm_mean_inferrer(gm, K, N, infer_method);
frmwork.add_func('update_u', mean_inferrer);

labeler = fmm_labeler( ...
    @(U, X) loglik(gm, U, X), ...
    struct('name', 'u', 'type', 'double', 'size', d), ...
    struct('name', 'x', 'type', 'double', 'size', d), ...
    K, N, infer_method);
frmwork.add_func('update_z', labeler);

frmwork.add_func('get_x', smi_store(X));

if has_u_pri
    frmwork.add_func('get_mu', smi_store(mu));
end

switch method
    case 'em'
        init_labeler = rand_labeler(K, N, 'b');
    case 'gs'
        init_labeler = rand_labeler(K, N, 'v');        
end

frmwork.add_func('init_z', init_labeler);

% add steps

frmwork.add_step('get_x', [], {'value', 'X'}, 1, 'init');
if has_u_pri
    frmwork.add_step('get_mu', [], {'value', 'mu'}, 1, 'init'); 
end
frmwork.add_step('init_z', [], {'result', 'Z'}, 1, 'init');    

inputs_for_u = {'x', 'X'; 'z', 'Z'};
if has_u_pri
    inputs_for_u = [inputs_for_u; {'mu', 'mu'}];
end
if est_Cx
    inputs_for_u = [inputs_for_u; {'Cx', 'Cx'}];
end
   
frmwork.add_step('update_u', inputs_for_u, {'u', 'U'});
    
frmwork.add_step('update_z', {'u', 'U'; 'x', 'X'}, {'Labels', 'Z'});    
 
% compile into a program

prg = smi_program.compile(frmwork);


