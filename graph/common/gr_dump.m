function gr_dump(G)
% Dumps a graph
%
%   gr_dump(G);
%

% Created by Dahua Lin, on Oct 31, 2010
%

%% verify input

if ~(isstruct(G) && isfield(G, 'tag'))
    error('gr_dump:invalidarg', 'G should be a struct representing a graph.');
end

if strcmp(G.tag, 'gr_edgelist')
    is_edgelist = true;
    is_adjlist = false;
    
elseif strcmp(G.tag, 'gr_adjlist')
    is_edgelist = true;
    is_adjlist = true;
    
else
    error('gr_dump:invalidarg', 'The input struct is invalid.');
end
   

%% main

fprintf('Basic Info:\n');
fprintf('------------------------\n');
fprintf('    # nodes = %d\n', G.n);
fprintf('    # edges = %d\n', G.m);
fprintf('    # direction = %s\n', direction_name(G));
fprintf('\n');

if is_edgelist
    fprintf('Edge List:\n');
    fprintf('----------------------\n');
    
    if isempty(G.w)
        for i = 1 : G.m
            fprintf('    [%d]: (%d, %d)\n', i, G.s(i)+1, G.t(i)+1);
        end
    else
        for i = 1 : G.m
            fprintf('    [%d]: (%d, %d) = %g\n', i, G.s(i)+1, G.t(i)+1, G.w(i));
        end
    end
    fprintf('\n');
end

if is_adjlist
    fprintf('Out neighbors:\n');
    fprintf('----------------------\n');
    for i = 1 : G.n
        b = G.o_os(i);
        d = G.o_ds(i);
        
        fprintf('    [%d] (deg = %d): ', i, d);
        
        for j = 1 : d
            e = G.o_es(b+j);
            
            if e < G.m
                fprintf('%d>%d ', e+1, G.o_ns(b+j)+1);
            else
                fprintf('%d<%d ', e-G.m+1, G.o_ns(b+j)+1);
            end
        end
        fprintf('\n');
    end
    fprintf('\n');
end

%% subfunctions

function s = direction_name(G)

if isfield(G, 'dty') && ~isempty(G.dty)       
    switch G.dty
        case 'd'
            s = 'directed';
        case 'u'
            s = 'undirected';
        otherwise
            s = 'INVALID!';
    end                    
else
    s = 'n/a';
end





