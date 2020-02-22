function [varargout] = model_basic(varargin)
% model_basic: mrna-protein constitutive gene expression, mass action, Simbiology. 
% This function runs a simple model and generates a plot, and optionally 
% returns a simbiology model, the corresponding exported model object, 
% an mcmc_info struct (with the ordered estimated species and parameters
% appended as a field too. 
% 
% INPUTS - Name value pairs, all optional. (defaults kick in if these are
% not specified)
% 1) 'mcmc_info', VALUE, where VALUE is a valid mcmc_info struct. See documentation
% in the file !TODO <write this documentation> for how to specify this. You
% can also see how it was specified within this fucntion (as a default). 
% 2) 'simtime', VALUE, where t_end is a numeric (double) specifying the time
% to the end of the simulation. 
% 
% OUTPUTS - generates a simbiology plot, and returns the following optional
% arguments in order:
% m: a model object, 
% em: an exported model object, 
% mi: an mcmc info struct. 
% 

%
% THIS MODEL:
% P + D <-> PD (kPDf, kPDr)
% PD -> P + D + m (kTX)
% m + R <-> mR (kmRf, kmRr)
% mR -> m + R + G (kTL)
% G -> null (kGd)
% m -> null (kmd)
% 
% 

% setup some defaults (defining the experiment structure suing the
% mcmc_info structure)
% setup the default mcmc_info object
estNames = {'kPDr'
    'kTX'
    'kmRr'
    'kTL'
    'kGd'
    'kmd'
    'P' 
    'R'};
paramranges = log([0.5 50
    0.005 5
    0.5 50
    0.005 5
    0.005 5
    0.005 5
    0.1 500
    0.1 500]);

dosedNames = {'D'}; %to be really sneaky, can try to dose mrna as well. 
% Fairly intractable experimentally though.  

dosedVals = [1]; % variations include multiple concentrations: 
% [0.5, 1], [0.1, 1, 10], [0.1, 1, 10, 100] and a study of how 
% identifiability changes with these variations. 

measuredSpecies = {{'m', 'mR'}, 'G'}; % the first species is a sum of the 
% ribosome bound RNA and free RNA, and the second is just the GFP. 

stdev = 1; % the ideal value of this is investigated in the project file 
% titled proj_basic_model.m

tightening = 1; % i have no idea what a good value is
nW = 300; % actual: 200 - 600 ish
stepsize = 1.5; % actual: 2 to 4 ish
niter = 30; % actual: 2 - 20 ish,
npoints = 2e4; % actual: 1e5 ish
thinning = 10; % actual: 10 to 40 ish

mi = struct(...
    'names_unord', {estNames}, ...
    'paramranges', {paramranges},...
    'dosednames', {dosedNames},...
    'dosedvals', {dosedVals},...
    'measuredspecies', {measuredSpecies}, ...
    'stdev', {stdev}, ...
    'tightening', {tightening}, ...
    'nW', {nW}, ...
    'stepsize', {stepsize}, ...
    'niter', {niter}, ...
    'npoints', {npoints}, ...
    'thinning', {thinning}, ...
    'parallel', true);

p = inputParser ;
p.addParameter(p, 'mcmc_info', mi);
p.addParameter(p, 'simtime', 200);

m1 = sbiomodel('constitutive_expression');

r1 = addreaction(m1,'P + D <-> PD');
Kobj = addkineticlaw(r1,'MassAction');
Kobj.ParameterVariableNames = {'kPDf','kPDr'};
addparameter(m1, 'kPDf', 1)
addparameter(m1, 'kPDr', 5)

r2 = addreaction(m1,'PD -> P + D + m');
Kobj = addkineticlaw(r2,'MassAction');
Kobj.ParameterVariableNames = {'kTX'};
addparameter(m1, 'kTX', 1)

r3 = addreaction(m1,'m + R <-> mR');
Kobj = addkineticlaw(r3,'MassAction');
Kobj.ParameterVariableNames = {'kmRf','kmRr'};
addparameter(m1, 'kmRf', 1)
addparameter(m1, 'kmRr', 5)

r4 = addreaction(m1,'mR -> m + R + G');
Kobj = addkineticlaw(r4,'MassAction');
Kobj.ParameterVariableNames = {'kTL'};
addparameter(m1, 'kTL', 1)

r5 = addreaction(m1,'G -> null');
Kobj = addkineticlaw(r5,'MassAction');
Kobj.ParameterVariableNames = {'kGd'};
addparameter(m1, 'kGd', 1)

r6 = addreaction(m1,'m -> null');
Kobj = addkineticlaw(r6,'MassAction');
Kobj.ParameterVariableNames = {'kmd'};
addparameter(m1, 'kmd', .1)


pol = sbioselect(m1, 'name', 'P');
pol.InitialAmount = 10;

dna = sbioselect(m1, 'name', 'D');
dna.InitialAmount = 5;

ribo = sbioselect(m1, 'name', 'R');
ribo.InitialAmount = 50;
m1.species
%%

cs = getconfigset(m1, 'active');
        set(cs, 'StopTime', p.simtime);
        
sd = sbiosimulate(m1);
sbioplot(sd);


%% DONT TOUCH THIS SECTION
% export and accelerate simbiology model object using estimated species
% and dosing species names

% select the parameters and species objects using the name array
ep = sbioselect(m1, 'Type', 'parameter', 'Name', ...
    mi.names_unord);% est parameters

es = sbioselect(m1, 'Type', 'species', 'Name', ...
    mi.names_unord);% est species

aps = [ep; es]; % active parameters and species

% reorder the parameter and species so they are in the same order as that
% in the model. 
eno = cell(length(aps), 1);% est names ordered

for i = 1:length(aps)
    eno{i} = aps(i).Name;
end

% 
ds = sbioselect(m1, 'Type', 'species', 'Name', mi.dosednames);

emo = export(m1, [ep; es; ds]); % exported model object, dosed species names. 
SI = emo.SimulationOptions;
SI.StopTime = p.simtime;
accelerate(emo);

mi.names_ord = eno;
mi.emo = emo; % !TODO: remove ?

% process variable outputs
nout = nargout;
outstuff = {m1, emo, mi};
varargout(1:nout) = outstuff(1:nout); % REMOVE: i think this should work. 

end

