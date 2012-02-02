function [w, b] = lsvm_dlm_train(S, hw, initsol, varargin)
%LSVM_DLM_TRAIN Linear SVM training via direct loss minimization
%
%   [w, b] = LSVM_DLM_TRAIN(S, hw, initsol, ...);
%
%       Learns a linear SVM model via direct minimization of an 
%       approximate regularized loss function, as
%
%           (1/2) ||w||^2 + sum_i c_i L( y_i * (w'*x_i+b) ).
%
%       Here, L is an approximate hinge loss, defined by
%
%           L(u) = 0                        when u >= 1 + h
%                = (1 + h - u)^2 / (4*h)    when 1 - h < u < 1 + h
%                = 1 - u                    when u <= 1 - h
%
%       Here, h is called the hinge transition width.
%
%       Input arguments:
%
%       - S:        The svm-problem struct 
%       - hw:       The hinge transition width [0 < hw < 1]
%       - initsol:  Initial guess of the solution, in form of [w; b]
%                   One can set initsol to [], then the function will
%                   be zeros(d+1, 1) as initial guess by default.
%
%       Output arguments:
%
%       - w:        The weight vector, size: d x 1
%       - b:        The offset scalar
%       
%       One can control the optimization procedure, by specifying the
%       following options in form of name/value pairs. See the help 
%       of newtonfmin for detailed information on what options are
%       available.
%       

% Created by Dahua Lin, on Jan 19, 2011
%

%% verify input arguments

if ~(is_svm_problem(S) && strcmp(S.type, 'class') && ...
        strcmp(S.kernel_type, 'linear'))
    error('lsvm_dlm_train:invalidarg', ...
        'S should be an svm-problem for classification using linear kernel.');
end

if ~(isfloat(hw) && isreal(hw) && hw > 0 && hw < 1)
    error('lsvm_dlm_train:invalidarg', ...
        'hw should be a real scalar in (0, 1).');
end

d = size(S.X, 1);
if nargin < 3 || isempty(initsol)
    s0 = zeros(d+1, 1);
else
    if ~(isfloat(initsol) && isreal(initsol) && isequal(size(initsol), [d 1]))
        error('lsvm_dlm_train:invalidarg', ...
            'initsol should be a real vector of size d x 1.');
    end
    s0 = initsol;
end

%% main

Xa = [S.X; ones(1, S.n)];
y = S.y;
c = S.C;
f = @(s) lsvm_dlm_objfun(s, Xa, y, c, hw);
sol = newtonfmin(f, s0, varargin{:});

w = sol(1:d);
b = sol(d+1);

%% objective function 

function [v, g, H] = lsvm_dlm_objfun(s, Xa, y, c, h)

d = size(Xa, 1) - 1;
w = s(1:d);
b = s(d+1);
b_reg = 1e-9;
c_sca = isscalar(c);

% evaluate linear predictor

u = y .* (s' * Xa);

% categorize samples

ti = find(abs(u - 1) < h);
wi = find(u <= 1 - h);

has_t = ~isempty(ti);
has_w = ~isempty(wi);

% objective value

if has_t
    ut = u(ti);
    rt = 1 + h - ut;
    if c_sca
        Lt = (1/(4*h)) * ( sum(rt.^2) * c );
    else
        ct = c(ti);    
        Lt = (1/(4*h)) * ( (rt.^2) * ct );
    end        
else
    Lt = 0;
end

if has_w
    uw = u(wi);
    
    if c_sca
        Lw = sum(1 - uw) * c;
    else
        cw = c(wi);
        Lw = (1 - uw) * cw;
    end
else
    Lw = 0;
end

v = 0.5 * ((w' * w) + b_reg * (b^2)) + (Lt + Lw);

% gradient

if nargout < 2; return; end

g = [w; b_reg * b];
if has_t
    gt = (1 / (2*h)) * rt;
    Xt = Xa(:, ti);
    if c_sca
        g = g - Xt * (c * (y(ti) .* gt).');
    else
        g = g - Xt * (ct .* (y(ti) .* gt).');
    end
end

if has_w
    Xw = Xa(:, wi);    
    if c_sca
        g = g - Xw * (c * y(wi)).';
    else
        g = g - Xw * (cw .* y(wi).');
    end
end

% Hessian

if nargout < 3; return; end

H = diag([ones(d, 1); b_reg]);
if has_t
    if c_sca
        H = H + (c/(2*h)) * (Xt * Xt');
    else
        H2 = Xt * bsxfun(@times, Xt', ct);
        H2 = 0.5 * (H2 + H2');
        H = H + (1/(2*h)) * H2;
    end
end


