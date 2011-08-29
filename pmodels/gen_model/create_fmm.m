function frmwork = create_fmm(model, obs, K, hparams, op)
% Creates a finite mixture model framework
%
%   frmwork = create_fmm(model, obs, K, hparams, op);
%       Creates a finite mixture model framework.
%
%       Input arguments:
%       - model:        the generative model (of class gen_model)
%       - obs:          the observations
%       - K:            the number of models to be in the mixture
%       - hparams:      the struct of hyper-parameters.
%                       Each field of the struct corresponds to a 
%                       hyper-parameter, and the field name is the
%                       name of hyper parameter specified in model.
%       - op:           The inference option: 
%                       - 'mapest': MAP estimation
%                       - 'sample': Gibbs sampling
%
%   Note: labeling prior is currently not incorporated yet.
%

% Created by Dahua Lin, on Aug 28, 2011
%

%% verify input arguments

if ~isa(model, 'gen_model')
    error('create_fmm:invalidarg', ...
        'The model should be of class gen_model or its derived class.');
end

prinfo = model.products_info;
painfo = model.params_info;
hpinfo = model.hyper_params_info;

m = numel(prinfo);  assert(m == 1);
np = numel(painfo); assert(np == 1);
nh = numel(hpinfo);

N = size(obs, ndims(obs));

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K > 0)
    error('create_fmm:invalidarg', ...
        'K should be a positive integer scalar.');
end

if ~(isstruct(hparams) && isscalar(hparams))
    error('create_fmm:invalidarg', 'hparams should be a struct scalar.');
end

if ~ischar(op)
    error('create_fmm:invalidarg', 'op should be a char string.');
end


%% main

% prepare functions

param_inferrer = fmm_param_inferrer(model, K, N, op);

labeler = fmm_labeler(model, K, N, op);


% build framework

frmwork = smi_frmwork();

% add vars

for i = 1 : m    
    frmwork.add_var(prinfo(i).name, prinfo(i).type, [prinfo(i).size, N]);
end

for i = 1 : np
    frmwork.add_var(painfo(i).name, painfo(i).type, [painfo(i).size, K]);
end

for i = 1 : nh
    frmwork.add_var(hpinfo(i).name, hpinfo(i).type, hpinfo(i).size);
end

label_slot = labeler.output_slots(1);
frmwork.add_var(label_slot.name, label_slot.type, label_slot.size);

loglik_slot = labeler.output_slots(2);
frmwork.add_var(loglik_slot.name, loglik_slot.type, loglik_slot.size);


% add functions

frmwork.add_func('update_param', param_inferrer);
frmwork.add_func('update_label', labeler);

frmwork.add_func([prinfo.name, '_store'], smi_store(obs));

hpa = false(1, nh);

for i = 1 : nh
    hpi = hpinfo(i);
    if isfield(hparams, hpi.name)
        hpv = hparams.(hpi.name);
        if ~isempty(hpv)
            fname = [hpi.name, '_store'];
            frmwork.add_func(fname, smi_store(hpv));
            hpa(i) = 1;
        end
    end
end

% add pre-steps

obs_step.func_name = [prinfo.name '_store'];
obs_step.slots = { 'value' };
obs_step.vars = { prinfo.name };

frmwork.add_init_steps(obs_step);

hyp_steps = cell(nh, 1);
for i = 1 : nh
    if hpa(i)
        s = [];
        s.func_name = [hpinfo(i).name '_store'];
        s.slots = { 'value' };
        s.vars = { hpinfo(i).name };
        hyp_steps{i} = s;
    end        
end
hyp_steps = [hyp_steps{hpa}];

frmwork.add_init_steps(hyp_steps);


% add steps

hpnames = {hpinfo.name};
hpnames = hpnames(hpa);
panames = {painfo.name};
prnames = {prinfo.name};

label_vname = label_slot.name;
loglik_vname = loglik_slot.name;

upa_step.func_name = 'update_param';
upa_step.slots = [hpnames, prnames, {label_vname}, panames];
upa_step.vars = upa_step.slots;

ula_step.func_name = 'update_label';
ula_step.slots = [panames, prnames, {label_vname}, {loglik_vname}];
ula_step.vars = ula_step.slots;

frmwork.add_steps( [ upa_step, ula_step ] );

% compile

frmwork.compile();


