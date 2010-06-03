function R = rotmat2(t)
% Gets the rotation matrix given radius values
%
%   R = rotmat2(t)
%       Returns the rotation matrix [cos(t) -sin(t); sin(t) cos(t)].
%       If t contains m elements, then R is a 2 x 2 x m array, with
%       R(:,;,i) corresponding to t(i).
%

% Created by Dahua Lin, on Apr 7, 2010
%

%% verify input arguments

if ~(isfloat(t) && isvector(t))
    error('rotmat2:invalidarg', 't should be a numeric vector.');
end

m = numel(t);

%% compute

if size(t, 1) > 1
    t = t.';
end
c = cos(t);
s = sin(t);

if m == 1
    R = [c -s; s c];
else
    R = reshape([c; s; -s; c], [2 2 m]);
end

