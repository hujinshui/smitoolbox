function write_mps(P, filename, fmt)
% Write an LP or QP optimization problem to MPS file
%
%   write_mps(P, filename);
%       writes an LP problem to the specified MPS file.
%
%   write_mps(P, filename, fmt);
%       writes an QP problem to the specified MPS file using specified
%       extension.
%
%       Note that the standard MPS format does not support QP in itself.
%       Different vendors have their respective extention to MPS format
%       to support QP and other problems.
%
%       The argument fmt here specifies which extension to use, whose
%       value can be
%       - 'cplex':  the ILOG CPLEX extensions.
%       - 'gurobi': the GUROBI extensions.
%       - 'mosek':  the MOSEK extensions.
%

%   Created by Dahua Lin, on April 16, 2011
%

%% verify input arguments

if ~(isstruct(P) && isfield(P, 'type') && ...
        (strcmp(P.type, 'lp') || strcmp(P.type, 'qp')) )
    error('write_mps:invalidarg', ...
        'P should be a lp_problem or qp_problem struct.');
end

is_qp = strcmp(P.type, 'qp');
if is_qp 
    if nargin < 3
        error('write_mps:invalidarg', ...
            'The 3rd argument (fmt) should be specified for QP problem.');
    end
    
    if ~ischar(fmt)
        error('write_mps:invalidarg', 'fmt should be a char string.');
    end
    
    FMT_CPLEX = 1;
    FMT_GUROBI = 2;
    FMT_MOSEK = 3;
    
    switch lower(fmt)
        case 'cplex'
            fmt_id = FMT_CPLEX;
        case 'gurobi'
            fmt_id = FMT_GUROBI;
        case 'mosek'
            fmt_id = FMT_MOSEK;
        otherwise
            error('write_mps:invalidarg', 'Unknown fmt name %s', fmt);
    end    
end


%% main

% open file

fid = fopen(filename, 'w');
if fid <= 0
    error('write_mps:ioerror', 'Failed to open file %s to write', filename);
end

% NAME section

[pastr, name] = fileparts(filename); %#ok<ASGLU>
fprintf(fid, '%-14s%s\n', 'NAME', upper(name));

% ROWS section

cons_names = write_rows_section(fid, P);

% COLUMNS section

var_names = write_columns_section(fid, P, cons_names);
fprintf(fid, '\n');


% RHS section

write_rhs_section(fid, P, cons_names);

% BOUNDS section

write_bounds_section(fid, P, var_names);

% Quadraic section (for QP)

if is_qp
    switch fmt_id
        case FMT_CPLEX
            write_quad_obj_cplex(fid, P.H, var_names);
        case FMT_GUROBI
            write_quad_obj_gurobi(fid, P.H, var_names);
        case FMT_MOSEK
            write_quad_obj_mosek(fid, P.H, var_names);
    end
end


% ENDATA

fprintf(fid, 'ENDATA\n\n');
    
% close file

fclose(fid);



%% standard sections

function cons_names = write_rows_section(fid, P)

fprintf(fid, 'ROWS\n');

fprintf(fid, ' N  COST\n');

ge_names = [];
le_names = [];
eq_names = [];

if ~isempty(P.A)          
    m1 = size(P.A, 1);
    ge_name_pat = sprintf('G%%0%dd', length(int2str(m1)));
    le_name_pat = sprintf('L%%0%dd', length(int2str(m1)));    
    ge_names = arrayfun(@(x) sprintf(ge_name_pat, x), (1:m1)', 'UniformOutput', false);
    le_names = arrayfun(@(x) sprintf(le_name_pat, x), (1:m1)', 'UniformOutput', false);
    
    for i = 1 : m1
        if is_val(P.bl, i)
            fprintf(fid, make_row('G', ge_names{i}));
        end
        if is_val(P.bu, i)
            fprintf(fid, make_row('L', le_names{i}));
        end
    end
end

if ~isempty(P.Aeq)
    m0 = size(P.Aeq, 1);
    eq_name_pat = sprintf('E%%%dd', length(int2str(m0)));    
    eq_names = arrayfun(@(x) sprintf(eq_name_pat, x), (1:m0)', 'UniformOutput', false);
    
    for i = 1 : m0
        fprintf(fid, make_row('E', eq_names{i}));
    end
end

cons_names.ge = ge_names;
cons_names.le = le_names;
cons_names.eq = eq_names;


function var_names = write_columns_section(fid, P, cons_names)

fprintf(fid, 'COLUMNS\n');

n = P.d;
cname_pat = sprintf('X%%0%dd', length(int2str(n)));
cnames = arrayfun(@(x) sprintf(cname_pat, x), (1:n).', 'UniformOutput', false);

ge_names = cons_names.ge;
le_names = cons_names.le;
eq_names = cons_names.eq;

for j = 1 : n
    
    cname = cnames{j};
    
    fv = 0;
    if ~isempty(P.f)
        fv = P.f(j);
    end
    
    if fv == 0    
        fprintf(fid, make_entry(cname));
    else
        fprintf(fid, make_entry(cname, 'COST', fv));
    end
    
    if ~isempty(P.A)
        [i, v] = find(P.A(:,j));
        
        if ~isempty(i)
            for k = 1 : numel(i)
                bl = is_val(P.bl, i(k));
                bu = is_val(P.bu, i(k));
                
                if bl && bu
                    fprintf(fid, make_entry(cname, ge_names{i(k)}, v(k), le_names{i(k)}, v(k)));
                elseif bl
                    fprintf(fid, make_entry(cname, ge_names{i(k)}, v(k)));
                elseif bu
                    fprintf(fid, make_entry(cname, le_names{i(k)}, v(k)));
                end                
            end
        end
    end
    
    if ~isempty(P.Aeq)
        [i, v] = find(P.Aeq(:,j));
        
        if ~isempty(i)
            for k = 1 : numel(i)
                fprintf(fid, make_entry(cname, eq_names{i(k)}, v(k)));
            end
        end
    end    
    
end

var_names = cnames;


function write_rhs_section(fid, P, cons_names)

fprintf(fid, 'RHS\n');

ge_names = cons_names.ge;
le_names = cons_names.le;
eq_names = cons_names.eq;

if ~isempty(P.A)
    m1 = size(P.A, 1);
    for i = 1 : m1;
        bl = is_val(P.bl, i);
        bu = is_val(P.bu, i);

        if bl && bu
            fprintf(fid, make_entry('RH', ge_names{i}, P.bl(i), le_names{i}, P.bu(i)));
        elseif bl
            fprintf(fid, make_entry('RH', ge_names{i}, P.bl(i)));
        elseif bu
            fprintf(fid, make_entry('RH', le_names{i}, P.bu(i)));
        end
    end
end

if ~isempty(P.Aeq)
    m0 = size(P.Aeq, 1);
    for i = 1 : m0
        fprintf(fid, make_entry('RH', eq_names{i}, P.beq(i)));
    end    
end


function write_bounds_section(fid, P, var_names)

fprintf(fid, 'BOUNDS\n');
n = P.d;

if isempty(P.l)
    if isempty(P.u)
        for j = 1 : n
            fprintf(fid, make_bnd('FR', var_names{j}));
        end
    else
        for j = 1 : n
            if is_val(P.u, j)
                fprintf(fid, make_bnd('UP', var_names{j}, P.u(j)));
            end
        end
    end
else
    if isempty(P.u)
        for j = 1 : n
            if is_val(P.l, j)
                fprintf(fid, make_bnd('LO', var_names{j}, P.l(j)));
            end
        end
    else
        for j = 1 : n
            if is_val(P.l, j)
                fprintf(fid, make_bnd('LO', var_names{j}, P.l(j)));
            end
            if is_val(P.u, j)
                fprintf(fid, make_bnd('UP', var_names{j}, P.u(j)));
            end
        end
    end
end

%% vendor-specific extended sections


function write_quad_obj_cplex(fid, H, var_names)

fprintf(fid, 'QMATRIX\n');

[I, J, V] = find(H);
si = find(I <= J);
I = I(si);
J = J(si);
V = V(si);

for k = 1 : numel(I)
    
    i = I(k);
    j = J(k);
    v = V(k);
    
    if i == j
        fprintf(fid, make_entry(var_names{i}, var_names{j}, v));
    else
        fprintf(fid, make_entry(var_names{i}, var_names{j}, v));
        fprintf(fid, make_entry(var_names{j}, var_names{i}, v));
    end    
end


function write_quad_obj_gurobi(fid, H, var_names)

fprintf(fid, 'QUADOBJ\n');

[I, J, V] = find(H);
si = find(I <= J);
I = I(si);
J = J(si);
V = V(si);

for k = 1 : numel(I)
    
    i = I(k);
    j = J(k);
    v = V(k);
    
    fprintf(fid, make_entry(var_names{i}, var_names{j}, v));   
end


function write_quad_obj_mosek(fid, H, var_names)

fprintf(fid, '%15s%s\n', 'QSECTION', 'COST');

[I, J, V] = find(H);
si = find(I <= J);
I = I(si);
J = J(si);
V = V(si);

for k = 1 : numel(I)
    
    i = I(k);
    j = J(k);
    v = V(k);
    
    if i == j
        fprintf(fid, make_entry(var_names{i}, var_names{j}, v));
    else
        fprintf(fid, make_entry(var_names{i}, var_names{j}, 2 * v));
    end
end


%% auxiliary formatting functions

function line = make_row(type, name)

line = sprintf(' %c  %s\n', type(1), name);


function line = make_entry(n0, n1, v1, n2, v2)

pb = '    ';

switch nargin
    case 1
        line = [pb, n0, '\n'];
    case 2
        line = [pb, sprintf('%-10s%s\n', n0, n1)];
    case 3
        line = [pb, sprintf('%-10s%-10s%g\n', n0, n1, v1)];
    case 5
        line = [pb, sprintf('%-10s%-10s%-15g%-10s%g\n', n0, n1, v1, n2, v2)];
end


function line = make_bnd(type, cname, v)

if nargin == 2
    line = sprintf(' %-2s %-10s%s\n', type, 'BND', cname);
else
    line = sprintf(' %-2s %-10s%-10s%g\n', type, 'BND', cname, v);
end

%% Other auxiliary function

function tf = is_val(a, i)

tf = ~isempty(a) && ~isinf(a(i)) && ~isnan(a(i));





