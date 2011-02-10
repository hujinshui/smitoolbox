function [vs, cs, info] = gatrwa(GM, ep, varargin)
% Solve tree reweighted approximation of Gaussian MRF
%
%   [vs, cs] = gatrwa(GM, ep, ...)
%   [vs, cs, info] = gatrwa(GM, ep, ...);
%   
%       solves the tree-reweighted approximation of Gaussian MRF given 
%       by GM.
%
%       Input arguments:
%       - GM:       an sgmrf object to represent a Gaussian Markov random
%                   field.
%       - ep:       the edge appearance probabilities. It can be either
%                   a scalar to assign the same value to all edges, or
%                   a vector of size m x 1. Here, m == GM.ne.
%
%       One can specify further options to control the solving process:
%
%       - 'sigma0':     initial guess of sigma (sqrt of variance)
%       - 'rho0':       initial guess of correlation coefficients
%       - 'maxiter':    the maximum number of iterations
%                       (default = 5000)
%       - 'tolfun':     the tolerance of objective change
%                       at convergence (default = 1e-8)
%       - 'tolx':       the tolerance of sigma change at
%                       convergence (default = 1e-6)
%       - 'order':      the order of updating
%                       (default = 1 : nv)
%       - 'display':    whether to display procedural info
%                       (default = false)
%
%       Output arguments:
%       - vs:           an nv x 1 vector that gives the marginal variance
%                       of each node.
%       - cs:           an ne x 1 vector that gives the marginal covariance
%                       between each pair of nodes connecting by edges.
%       - info:         an information struct that contains the following
%                       fields about the procedure:
%
%                       - 'niters':     the number of elapsed iterations
%                       - 'converged':  whether the procedure converged
%                       - 'energy':     the expected energy 
%                                       ( = (1/2) * tr(J * C) )
%                       - 'entropy':    the approximated entropy
%                       - 'objv':       objective value: energy - entropy
%

%   History
%   -------
%       - Created by Dahua Lin, on Feb 9, 2011
%


%% verify input arguments

if ~isa(GM, 'sgmrf')
    error('gatrwa:invalidarg', 'GM should be an sgmrf object.');
end

if ~(isfloat(ep) && (isscalar(ep) || isequal(size(ep), [GM.ne, 1])))
    error('gatrwa:invalidarg', 'ep should be an ne x 1 numeric vector.');
end
if isscalar(ep)
    ep = constmat(GM.ne, 1, ep);
end

opts = struct( ...
    'sigma0', [], ...
    'rho0', [], ...
    'maxiter', 5000, ...
    'tolfun', 1e-8, ...
    'tolx', 1e-6, ...
    'order', 1:GM.nv, ...
    'display', false);

if ~isempty(varargin)
    opts = parse_options(opts, GM, varargin);
end


%% main

% initialize

if isempty(opts.sigma0)
    s = 1 ./ GM.Jdv;
else
    s = opts.sigma0;
end

if isempty(opts.rho0)
    r = zeros(GM.ne, 1);
else
    r = opts.rho0;
end

% main loop

it = 0;
converged = false;

energy = GM.qenergy_sr(s, r);
entropy = eval_entropy(ep, s, r);
objv = energy - entropy;

if opts.display
    fprintf('  Iter     obj.value     objv.ch    sigma.ch\n');
    fprintf('------------------------------------------------\n');
end

while ~converged && it < opts.maxiter
    it = it + 1;
    pre_objv = objv;
    pre_sig = s;
    
    [s, r] = cupdate(GM, ep, s, r, opts.order);
    
    energy = GM.qenergy_sr(s, r);
    entropy = eval_entropy(ep, s, r);
    objv = energy - entropy;
    
    chf = objv - pre_objv;
    chx = max(abs(s - pre_sig));
    converged = abs(chf) < opts.tolfun && chx < opts.tolx;
    
    if opts.display
        fprintf(' %5d  %12.4f  %10.3g  %10.3g\n', ...
            it, objv, chf, chx);
    end
end

if opts.display
    if converged
        fprintf('converged with %d iterations.\n', it);
    else
        fprintf('terminated without convergence.\n');
    end
end

% output

vs = s.^2;
cs = r .* (s(GM.es) .* s(GM.et));

if nargout >= 1
    info = struct( ...
        'niters', it, ...
        'converged', converged, ...
        'energy', energy, ...
        'entropy', entropy, ...
        'objv', objv);
end


%% Sub functions

function [s, r] = cupdate(GM, ep, s, r, ord)
% Perform combined updates
%
%   obj.cupdate(ord);
%
%       updates the sigma values of specified vertices along
%       with the rho values of their incident edges
%


ord = int32(ord - 1);
[s, r] = gatrwa_cimp(3, GM.Jgr, GM.Jdv, s, r, ep, ord);
        

function fv = eval_entropy(ep, s, r)

n = length(s);

hc = (log(2*pi) + 1) / 2;
H = sum(log(s)) + hc * n;
I = (-0.5) * (ep' * log(1 - r.^2));

fv = H - I;


%% Option parsing

function opts = parse_options(opts, GM, oplist)

onames = oplist(1:2:end);
ovals = oplist(2:2:end);

if ~(length(onames) == length(ovals) && iscellstr(onames))
    error('gatrwa:invalidarg', ...
        'The option list is invalid.');
end

for i = 1 : length(onames)
    cn = onames{i};
    cv = ovals{i};
    
    switch lower(cn)
        case 'sigma0'
            if ~(isfloat(cv) && isequal(size(cv), [GM.nv, 1]))
                error('gatrwa:invalidarg', ...
                    'sigma0 should be a numeric vector of size nv x 1.');
            end
            if ~isa(cv, 'double'); cv = double(cv); end
            opts.sigma0 = cv;
            
        case 'rho0'
            if ~(isfloat(cv) && isequal(size(cv), [GM.ne, 1]))
                error('gatrwa:invalidarg', ...
                    'rho0 should be a numeric vector of size ne x 1.');
            end
            if ~isa(cv, 'double'); cv = double(cv); end
            opts.rho0 = cv;
            
        case 'maxiter'
            if ~(isnumeric(cv) && isscalar(cv) && cv >= 1)
                error('gatrwa:invalidarg', ...
                    'maxiter should be a positive integer.');
            end
            opts.maxiter = cv;
            
        case 'tolfun'
            if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0)
                error('gatrwa:invalidarg', ...
                    'tolfun should be a positive real scalar.');
            end
            opts.tolfun = cv;
            
        case 'tolx'
            if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0)
                error('gatrwa:invalidarg', ...
                    'tolx should be a positive real scalar.');
            end
            opts.tolx = cv;
            
        case 'order'
            if ~(isnumeric(cv) && isvector(cv))
                error('gatrwa:invalidarg', ...
                    'order should be a numeric vector.');
            end
            opts.order = cv;
            
        case 'display'
            if ~(islogical(cv) && isscalar(cv))
                error('gatrwa:invalidarg', ...
                    'display should be a logical scalar.');
            end
            opts.display = cv;
            
        otherwise
            error('gatrwa:invalidarg', ...
                'Invalid option name %s', cn);
    end
end

