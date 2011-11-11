function [x, y, thres] = roc_curve(scores, L0, xname, yname, varargin)
% Computes the ROC curve from prediction results
%
%   [x, y, thres] = roc_curve(scores, L0, xname, yname, ...);
%
%       Computes the ROC curve from scores. 
%
%       Input arguments:
%       - scores:       The prediction scores.
%       - L0:           The ground-truth label of the scores. [logical]
%       - xname:        The name of quantity for x-axis
%       - yname:        The name of quantity for y-axis
%
%       Output arguments:
%       - x:            The x-coordinates of the sampled points at curve
%       - y:            The y-coordinates of the sampled points at curve
%       - thres:        The corresponding thresholds
%
%       Here, xname and yname can be chosen from either of the following:
%       - 'tp':         true-positive ratio = #[1->1] / #[*->1]
%       - 'fp':         false-positive ratio = #[0->1] / #[*->1]
%       - 'tn':         true-negative ratio = #[0->0] / #[*->0]
%       - 'fn':         false-negative ratio = #[1->0] / #[*->0]
%       - 'precision':  the retrieval precision = #[1->1] / #[*->1]
%       - 'recall':     the retrieval recall rate = #[1->1] / #[1->*]
%
%       Here, #[i->j] means the number of samples whose ground-truth
%       label is i and predicted label is j. i and j here can be either
%       0 (negative), 1 (positive), and * (both).
%
%       Procedure description:
%       - A sequence of thresholds are set up, with each threshold,
%         the score above it is considered as positive, and the remaining
%         are considered as negative. This procedure is implemented in
%         an efficient way through a multi-bin histogram of scores.
%       - According to these results, we compute the quantity specified
%         by xname and yname respectively.
%       - Via (linear) interpolation, we produce the y-values for the
%         given x-values.
%         
%       One can specify following options to customize the behavior:
%       - 'pred':       How to make the prediction, which can be
%                       'gt':   predict 1 when score > threshold; (default)
%                       'ge':   predict 1 when score >= threshold;
%                       'lt':   predict 1 when score < threshold;
%                       'le':   predict 1 when score <= threshold;
%
%       - 'sdensity':   The sample density: the number of thresholds
%                       to sample. (default = 100).
%                       It can also be a sorted vector that explicitly
%                       gives the sequence of thresholds.
%
%                       Note the number of thresholds must be at least 3.
%

%   History
%   -------
%       - Created by Dahua Lin, on April 24, 2011
%

%% verify input arguments

% basic args

if ~(isnumeric(scores) && isreal(scores))
    error('roc_curve:invalidarg', 'scores should be a numeric array.');
end

if ~(islogical(L0) && isequal(size(L0), size(scores)))
    error('roc_curve:invalidarg', ...
        'L0 should be a logical array of the same size as scores.');
end

if ~(ischar(xname) && ischar(yname))
    error('roc_curve:invalidarg', ...
        'xname and yname should be both strings.');
end
xfunc = name_to_cfunc(xname);
yfunc = name_to_cfunc(yname);


% options

opts.pred = 'gt';
opts.sdensity = 100;

if ~isempty(varargin)
    
    onames = varargin(1:2:end);
    ovals = varargin(2:2:end);
    
    if ~(numel(onames) == numel(ovals) && iscellstr(onames))
        error('roc_curve:invalidarg', 'The options are not correctly given.');
    end
    
    for i = 1 : numel(onames)
        cn = onames{i};
        cv = ovals{i};
        
        switch lower(cn)
            case 'pred'
                if ~(ischar(cv) && ismember(cv, {'gt', 'ge', 'lt', 'le'}))
                    error('roc_curve:invalidarg', ...
                        'The option pred is invalid.');
                end
            case 'interp'
                if ~(ischar(cv) && ...
                        ismember(cv, {'cubic', 'spline', 'linear', 'nearest'}))
                    error('roc_curve:invalidarg', ...
                        'The option interp is invalid.');
                end
            case 'sdensity'
                if (isnumeric(cv) && isreal(cv) && ...
                        (isscalar(cv) || isvector(cv)))
                    if isscalar(cv)
                        if cv < 3
                            error('roc_curve:invalidarg', ...
                                'sdensity must be at least 3.');
                        end
                    end
                else
                    error('roc_curve:invalidarg', ...
                        'The option sdensity is invalid.');
                end
            otherwise
                error('roc_curve:invalidarg', ...
                    'The option name %s is unknown', cn);            
        end
        
        opts.(cn) = cv;        
    end        
end

%% main

% sort and categorize score values

if ~(isvector(scores) && size(scores, 2) == 1)
    scores = scores(:);
    L0 = L0(:);
end
[scores, sord] = sort(scores);
L0 = L0(sord);

scores_t = scores(L0);
scores_f = scores(~L0);

clear sord;

% delimit scores into threshold bins

if isscalar(opts.sdensity)
    thres = linspace(scores(1), scores(end), opts.sdensity);
else
    thres = opts.sdensity;
end

switch opts.pred
    case 'gt'
        hdir = 'R';
        cdir = 1;
    case 'ge'
        hdir = 'L';
        cdir = 1;
    case 'lt'
        hdir = 'L';
        cdir = -1;
    case 'le'
        hdir = 'R';
        cdir = -1;
end
        
H_t = intcount([0, numel(thres)], pieces(scores_t, thres, hdir, 'sorted'));
H_f = intcount([0, numel(thres)], pieces(scores_f, thres, hdir, 'sorted'));


% calculate desired quantities

x = xfunc(H_t, H_f, cdir);
y = yfunc(H_t, H_f, cdir);

to_retain = ~isnan(x) & ~isnan(y);
x = x(to_retain);
y = y(to_retain);
thres = thres(to_retain);


%% Calculation functions

function f = name_to_cfunc(name)

switch lower(name)
    case 'tp'
        f = @calc_tp;
    case 'tn'
        f = @calc_tn;
    case 'fp'
        f = @calc_fp;
    case 'fn'
        f = @calc_fn;
    case 'precision'
        f = @calc_precision;
    case 'recall'
        f = @calc_recall;
    otherwise
        error('roc_curve:invalidarg', ...
            'The quantity name %s is invalid', name);
end

        
function v = calc_tp(H_t, H_f, cdir)

v = csum(H_t, cdir) ./ csum(H_t + H_f, cdir);


function v = calc_tn(H_t, H_f, cdir)

v = csum(H_f, -cdir) ./ csum(H_t + H_f, -cdir);


function v = calc_fp(H_t, H_f, cdir)

v = csum(H_f, cdir) ./ csum(H_t + H_f, cdir);


function v = calc_fn(H_t, H_f, cdir)

v = csum(H_t, -cdir) ./ csum(H_t + H_f, -cdir);


function v = calc_precision(H_t, H_f, cdir)

v = csum(H_t, cdir) ./ csum(H_t + H_f, cdir);


function v = calc_recall(H_t, H_f, cdir) %#ok<INUSL>

v = csum(H_t, cdir) / sum(H_t);



function cs = csum(h, dir)

if dir > 0
    cs = cumsum(h(end:-1:1));
    cs = cs(end-1:-1:1);
else
    cs = cumsum(h);
    cs = cs(1:end-1);
end


