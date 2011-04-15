function [x, fval, flag, info] = mosek_lpqp(P, options)
% The wrapper of MOSEK optimizer for solving LP or QP problem.
%
%   x = mosek_lpqp(P);
%   x = mosek_lpqp(P, options);
%
%   [x, fval]             = mosek_lpqp( ... );
%   [x, fval, flag]       = mosek_lpqp( ... );
%   [x, fval, flag, info] = mosek_lpqp( ... );
%
%   Use MOSEK LP/QP optimizer to solve the LP problem specified by the
%   lp_problem struct P.
%
%   options is a struct with MOSEK parameter fields. In addition, one 
%   can specify the following option:
%   - 'Display':    whose values can be 'off', 'final', or 'iter'.
%
%   Output arguments:
%
%       - x:    the solution
%       - fval: the value of objective function at x
%       - flag: the return code that describes the exit condition.
%       - info: the struct that contains information of the optimization
%               procedure and result.
%

% Created by Dahua Lin, on April 14, 2011
%

%% verify input arguments

if ~(isstruct(P) && isfield(P, 'type') && (strcmp(P.type, 'lp') || strcmp(P.type, 'qp')))
    error('mosek_lp:invalidarg', 'P should be a lp_problem or qp_problem struct.');
end

is_qp = strcmp(P.type, 'qp');


if nargin < 2 || isempty(options)
    options = [];
    options.Display = 'off';
else
    if ~isstruct(options)
        error('mosek_lpqp:invalidarg', 'options should be a MOSEK parameter struct.');
    end
end

[cmd, verbosity, param] = msksetup(true, options);


%% main

% prepare 

A = P.A;
bl = P.bl;
bu = P.bu;

Ae = P.Aeq;
be = P.beq;

a = [];
blc = [];
buc = [];
    
if ~isempty(A)
    if ~isempty(Ae)  % A & Ae
        a = [A; Ae];
        if isempty(bl)
            blc = [-inf(size(A,1), 1); be];
        else
            blc = [bl; be];
        end
        if isempty(bu)
            buc = [inf(size(A,1), 1); be];
        else
            buc = [bu; be];
        end
        
    else  % only A
        a = A;
        blc = bl;
        buc = bu;
        
    end
    
elseif ~isempty(Ae)  % only Ae
    a = Ae;
    blc = be;
    buc = be;
            
end

if ~issparse(a)
    a = sparse(a);
end


% setup problem

prob = [];
if is_qp
    [prob.qosubi,prob.qosubj,prob.qoval] = find(tril(sparse(P.H))); 
end

prob.c = P.f;
prob.a = a;
prob.blc = blc;
prob.buc = buc;
prob.blx = P.l;
prob.bux = P.u;

% solve

[rcode, res] = mosekopt(cmd, prob, param);

mskstatus('mosek_lpqp', verbosity, 0, rcode, res);

% extract output

if isfield(res, 'sol')
    x = res.sol.itr.xx;
else
    x = [];
end

if nargout >= 2
    if ~isempty(x)
        if ~is_qp
            fval = P.f' * x;
        else
            fval = (x' * P.H * x) / 2 + P.f' * x;
        end
    else
        fval = nan;
    end
end

if nargout >= 3
    flag = rcode;
end

if nargout >= 4
    info = res;
end


