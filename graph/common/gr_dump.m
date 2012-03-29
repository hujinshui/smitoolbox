function gr_dump(g)
% Dump the information of a graph
%
%   gr_dump(g);
%

% Created by Dahua Lin, on Oct 28, 2011
%

%% verify input

if ~is_gr(g)
    error('gr_dump:invalidarg', 'The input argument is not a valid graph struct.');
end

%% main

if g.dty == 'd'
    dtname = 'Directed';
else
    dtname = 'Undirected';
end

% basic information

fprintf('\n');
fprintf('%s graph with %d nodes and %d edges\n', dtname, g.n, g.m);
fprintf('----------------------------------------------\n');

% edges

fprintf('\n');
fprintf('Edges:\n');
if g.dty == 'd'
    for i = 1 : g.m
        fprintf('    %d --> %d\n', g.edges(1,i), g.edges(2,i));
    end
else
    for i = 1 : g.m
        fprintf('    %d --- %d\n', g.edges(1,i), g.edges(2,i));
    end
end

% neighborhood

if g.has_nbs    
    fprintf('\n');
    fprintf('Neighborhoods: \n');
    
    for i = 1 : g.n  
        deg = g.o_degs(i);
        co = g.o_os(i);
        
        fprintf('    [%d] (deg = %d): ', i, deg);
        for j = 1 : deg
            fprintf('%d ', g.o_nbs(co+j));
        end        
        fprintf('\n');
    end    
end
    
fprintf('\n');


