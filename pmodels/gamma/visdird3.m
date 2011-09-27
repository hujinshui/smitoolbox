function h = visdird3(alpha, form)
% Visualize a Dirichlet distribution with K == 3
%
%   visdird3(alpha);
%   visdird3(alpha, 'image');
%   visdird3(alpha, 'mesh');
%
%       visualizes a Dirichlet distribution with K = 3 as a colored 
%       image or mesh. (image by default);
%
%       alpha is either a scalar or a vector of length 3.
%
%   h = visdird3(alpha);
%   h = visdird3(alpha, 'image');
%   h = visdird3(alpha, 'mesh');
%
%       This statement returns the handle to the figure.
%

% Created by Dahua Lin, on Sep 17, 2011
%

%% verify inputs

if ~(isfloat(alpha) && isreal(alpha) && ...
        (isscalar(alpha) || (isvector(alpha) && numel(alpha) == 3)))
    error('visdird3:invalidarg', ...
        'alpha should be either a real scalar or a real vector of length 3.');
end

if size(alpha, 2) > 1
    alpha = alpha.';
end

if nargin < 2
    form = 'image';
else
    if ~(ischar(form) && (strcmpi(form, 'image') || strcmpi(form, 'mesh')))
        error('visdird3:invalidarg', ...
            'The 2nd argument should be either ''image'' or ''mesh''.');
    end
    form = lower(form);
end


%% main

% construct distribution

D = dirichletd(3, alpha, 'pre');

% gather the points at which pdf is calculated

tarref = [-1 0 1; 0 sqrt(3) 1; 1 0 1]';

switch form
    case 'image'
        nx = 480;
        ny = 300;
    case 'mesh'
        nx = 96;
        ny = 60;
end

tx = (1 : nx) / (nx + 1) * 2 - 1;
ty = (1 : ny) / (ny + 1) * sqrt(3);

[xx, yy] = meshgrid(tx, ty);

tar = [xx(:) yy(:)]';
np = size(tar, 2);

s = tarref \ [tar; ones(1, np)];
s = bsxfun(@times, s, 1 ./ sum(s, 1));

is_valid = all(s > 1e-3, 1);

% compute pdf

v = zeros(1, np);
v(is_valid) = D.pdf(s(:, is_valid));
vv = reshape(v, size(xx));

% visualize

switch form
    case 'image'
        h = imagesc(tx, ty', vv, [0, max(v)]);        
        axis xy;        
        axis([-1 1 0 sqrt(3)]);
        axis equal;        
        colorbar;
        set(gca, 'XTick', [], 'YTick', []);
        
    case 'mesh'
        h = mesh(xx, yy, vv);
        axis([-1 1 0 sqrt(3) 0 max(v)]);
        set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);       
end


