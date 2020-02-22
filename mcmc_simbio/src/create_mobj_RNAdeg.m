function [ Mobj1 ] = create_mobj_RNAdeg(extract, varargin)
% commit test
%create_mobj_geneexpr Create a standard TXTL toolbox gene expression
%circuit
%  [ Mobj1 ] = create_mobj_geneexpr(extract, speciesGroups, globalKdRules)
%
% Input Arguments: 
%
% extract = a string specifying the extract batch for the txtl Modeling
% toolbox, example: 'E30VNPRL'
%
% speciesGroups = a structure containing info on additional species to be
% created which are the sum of some of the existing species. The main
% example of this is TotalRNA, which is the sum of all the GFP RNA in
% the system. 
% 
% globalKdRules is a struct which looks like this:
% globalKdRules = 
% 
% 2x1 struct array with fields:
% 
%     rxStr
%     paramName
%     kdVal
%     fVal
% 
% globalKdRules(1)
%         rxStr: '[RNA rbs--deGFP] + Ribo <-> [Ribo:RNA...'
%     paramName: 'TXTL_RBS_Ribo_deGFP'
%         kdVal: 452.1739
%          fVal: 0.2300
% 
% See the file function_design_script for an example of the construction of
% the struct. 
% 
% OUTPUT Arguments
% Mobj1: a Simbiology model Object containing the constitutive gene
% expression circuit. 
%
if nargin > 1
    tspan = varargin{1};
    customtime = true;
end

% Create the standard model object
tube1 = txtl_extract(extract);
tube2 = txtl_buffer(extract);
tube3 = txtl_newtube('gene_expression');
txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', 0, 'plasmid');					
Mobj1 = txtl_combine([tube1, tube2, tube3]);
cs1 = getconfigset(Mobj1);
set(cs1.RuntimeOptions, 'StatesToLog', 'all');

set(cs1.SolverOptions, 'AbsoluteToleranceScaling', 1);
set(cs1.SolverOptions, 'AbsoluteTolerance', 1.0e-6);
set(cs1.SolverOptions, 'AbsoluteToleranceStepSize', tspan(end)*1.0e-6*0.1);
set(cs1.SolverOptions, 'RelativeTolerance', 1.0e-6);
% try: AbsoluteToleranceStepSize = AbsoluteTolerance * StopTime * 0.1
tic
if customtime
    [simData1] = txtl_runsim(Mobj1,tspan(end));
    else
[simData1] = txtl_runsim(Mobj1,8*60*60);
end
toc

%% Globalize irreverisble reactions too
rxstr = {'[protein deGFP] -> [protein deGFP*]'
    '[term_RNAP70:DNA p70--rbs--deGFP] -> RNAP70 + [DNA p70--rbs--deGFP]'
    '[RNA rbs--deGFP:RNase] -> RNase'
    '[AA:AGTP:Ribo:RNA rbs--deGFP:RNase] -> AGTP + AA + Ribo + RNase'
    '[Ribo:RNA rbs--deGFP:RNase] -> Ribo + RNase'};

pname_cat = {'TXTL_PROT_deGFP_MATURATION'
'TXTL_RNAPBOUND_TERMINATION_RATE'
'TXTL_RNAdeg_catalysis'
'TXTL_RNAdeg_catalysis'
'TXTL_RNAdeg_catalysis'};

k_cat = {0.00231049 
0.05
0.00277778
0.00277778
0.00277778};

globalCatalysisParams = struct('rxStr', rxstr,...
    'paramName', pname_cat,...
    'paramValue', k_cat);


for k = 1: length(globalCatalysisParams)
createFrateRules(Mobj1, globalCatalysisParams(k).rxStr,...
    globalCatalysisParams(k).paramName,...
    globalCatalysisParams(k).paramValue);
end

%% Define the parameters to move to global scope to allow for Kd estimation
rxstr = {'[RNA rbs--deGFP] + Ribo <-> [Ribo:RNA rbs--deGFP]';
    '[DNA p70--rbs--deGFP] + RNAP70 <-> [RNAP70:DNA p70--rbs--deGFP]'
    'RNAP + [protein sigma70] <-> RNAP70'
    '[RNAP70:DNA p70--rbs--deGFP] + AGTP <-> [AGTP:RNAP70:DNA p70--rbs--deGFP]'
    '[RNAP70:DNA p70--rbs--deGFP] + CUTP <-> [CUTP:RNAP70:DNA p70--rbs--deGFP]'
    '[AGTP:RNAP70:DNA p70--rbs--deGFP] + CUTP <-> [CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP]'
    '[CUTP:RNAP70:DNA p70--rbs--deGFP] + AGTP <-> [CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP]'
    '[Ribo:RNA rbs--deGFP] + AA + AGTP <-> [AA:AGTP:Ribo:RNA rbs--deGFP]'
    '[RNA rbs--deGFP] + RNase <-> [RNA rbs--deGFP:RNase]'
    '[AA:AGTP:Ribo:RNA rbs--deGFP] + RNase <-> [AA:AGTP:Ribo:RNA rbs--deGFP:RNase]'
    '[Ribo:RNA rbs--deGFP] + RNase <-> [Ribo:RNA rbs--deGFP:RNase]'};

pname = {'TXTL_RBS_Ribo_deGFP';
    'TXTL_P70_RNAPbound_deGFP'
    'TXTL_RNAP_S70'
    'TXTL_NTP_RNAP_1'
    'TXTL_NTP_RNAP_1'
    'TXTL_NTP_RNAP_2'
    'TXTL_NTP_RNAP_2'
    'TXTL_AA'
    'TXTL_RNAdeg'
    'TXTL_RNAdeg'
    'TXTL_RNAdeg'};
% using defaults from above. 
Kdval = {104/0.23;
    725.07/0.06;
    0.1/10;
    1.2e+10 / 100000;
    1.2e+10 / 100000;
    1.2e7 / 100;
    1.2e7 / 100;
    325046 / 9.055;
    2000 / 10;
    2000 / 10;
    2000 / 10}; 
Fval = {0.23;
    0.06;
    10;
    100000;
    100000;
    100;
    100;
    9.055;
    10;
    10;
    10}; 

globalKdRules = struct('rxStr', rxstr,...
    'paramName', pname,...
    'kdVal', Kdval,...
    'fVal', Fval);


% Create Kd rules, adding parameters to the global scope
for k = 1: length(globalKdRules)
createKdRules(Mobj1, globalKdRules(k).rxStr,...
    globalKdRules(k).paramName,...
    globalKdRules(k).kdVal,...
    globalKdRules(k).fVal);
end



%%
speciesGroups = struct('summedSpeciesName', {'TotalRNA'},...
    'speciesToSum', {{'[RNA rbs--deGFP]', '[Ribo:RNA rbs--deGFP]',...
    '[AA:AGTP:Ribo:RNA rbs--deGFP]', '[RNA rbs--deGFP:RNase]',...
    '[AA:AGTP:Ribo:RNA rbs--deGFP:RNase]', '[Ribo:RNA rbs--deGFP:RNase]'}});

% Create the rule for total RNA (or total whatever)
for i = 1:length(speciesGroups)
    expr = [speciesGroups(i).summedSpeciesName ' = '];
    nSpToSum = length(speciesGroups(i).speciesToSum);
    for j = 1:nSpToSum-1
        expr = [expr speciesGroups(i).speciesToSum{j} ' + '];
    end
    expr = [expr speciesGroups(i).speciesToSum{nSpToSum}]; 
    txtl_addspecies(Mobj1, speciesGroups(i).summedSpeciesName, 0);
    ruleObj = addrule(Mobj1, expr, 'repeatedAssignment');
end

%% Create a rule for consumption reaction using length based formula

%%%% FOR TX
% calculate tx rate. the rna len can be extracted from the userdata if
% desired. dont need to do that just yet. 
rnalen = 20 + 1000; %utr length and gene length in base pairs. 
elongrate = 1.5 ; % the parameter name is k_elon. this gets varied in estimation. 
rx = sbioselect(Mobj1, 'reaction', ...
    '[CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP] -> [term_RNAP70:DNA p70--rbs--deGFP] + [RNA rbs--deGFP]'); %
pTarget = sbioselect(rx, 'type', 'parameter', 'Name', 'TXTL_transcription_rate1');
pTarget.delete;
rx.KineticLaw.ParameterVariableNames = 'TXTL_transcription_rate1'; 

if isempty(sbioselect(Mobj1,'Type','Parameter', 'Name', 'TXTL_transcription_rate1'))
addparameter(Mobj1, 'TXTL_transcription_rate1', elongrate/rnalen);
end

if isempty(sbioselect(Mobj1,'Type','Parameter', 'Name', 'k_elon'))
addparameter(Mobj1, 'k_elon', elongrate);
end

ruleStr = ['TXTL_transcription_rate1 = k_elon/' num2str(rnalen)];

if isempty(sbioselect(Mobj1,'Type','Rule', 'Rule', ruleStr))
	addrule(Mobj1, ruleStr, 'initialAssignment');
end

rx = sbioselect(Mobj1, 'reaction', ...
    '[CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP] -> RNAP70 + [DNA p70--rbs--deGFP]'); %
pTarget = sbioselect(rx, 'type', 'parameter', 'Name', 'TXTL_NTP_consumption');
pTarget.delete;
rx.KineticLaw.ParameterVariableNames = 'TXTL_k_con_TX';

if isempty(sbioselect(Mobj1,'Type','Parameter', 'Name', 'TXTL_k_con_TX'))
addparameter(Mobj1, 'TXTL_k_con_TX', (rnalen/4-1)*(elongrate/rnalen));
end

ruleStr = ['TXTL_k_con_TX = (' num2str(rnalen) '/4-1)*(k_elon/' num2str(rnalen) ')'];
if isempty(sbioselect(Mobj1,'Type','Rule', 'Rule', ruleStr))
	addrule(Mobj1, ruleStr, 'initialAssignment');
end

%%%% FOR TL
genelen = 1000; %utr length and gene length in base pairs. 
elongrate_prot = 4 ; % the parameter name is k_elon. this gets varied in estimation. 

rx = sbioselect(Mobj1, 'reaction', ...
    '[AA:AGTP:Ribo:RNA rbs--deGFP] -> [RNA rbs--deGFP] + [protein deGFP] + Ribo'); %

pTarget = sbioselect(rx, 'type', 'parameter', 'Name', 'TXTL_TL_rate');
pTarget.delete;
rx.KineticLaw.ParameterVariableNames = 'TXTL_TL_rate'; 

if isempty(sbioselect(Mobj1,'Type','Parameter', 'Name', 'TXTL_TL_rate'))
addparameter(Mobj1, 'TXTL_TL_rate', elongrate_prot/genelen);
end

if isempty(sbioselect(Mobj1,'Type','Parameter', 'Name', 'k_elon_prot'))
addparameter(Mobj1, 'k_elon_prot', elongrate_prot);
end

ruleStr = ['TXTL_TL_rate = k_elon_prot/' num2str(genelen)];

if isempty(sbioselect(Mobj1,'Type','Rule', 'Rule', ruleStr))
	addrule(Mobj1, ruleStr, 'initialAssignment');
end

rx = sbioselect(Mobj1, 'reaction', ...
    '[AA:AGTP:Ribo:RNA rbs--deGFP] -> [RNA rbs--deGFP] + Ribo'); %

pTarget = sbioselect(rx, 'type', 'parameter', 'Name', 'TXTL_TL_AA_consumption');
pTarget.delete;
rx.KineticLaw.ParameterVariableNames = 'TXTL_k_con_TL'; 


if isempty(sbioselect(Mobj1,'Type','Parameter', 'Name', 'TXTL_k_con_TL'))
addparameter(Mobj1, 'TXTL_k_con_TL', (genelen-1)*(elongrate_prot/genelen));
end

ruleStr = ['TXTL_k_con_TL = (' num2str(genelen) '-1)*(k_elon_prot/' num2str(genelen) ')'];
if isempty(sbioselect(Mobj1,'Type','Rule', 'Rule', ruleStr))
	addrule(Mobj1, ruleStr, 'initialAssignment');
end

end

