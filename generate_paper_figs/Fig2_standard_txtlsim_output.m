% generate the constitutive gene expression figure, as described at the end
% of the tile analysis_vnprl_F2.m. 
clc

projdir = [pwd '/mcmc_simbio/projects/proj_vnprl/'];

ts = '20190223_024333_1_2';
tstamp = {ts};
nIterID = {20};
load([projdir '/simdata_' ts '/full_variable_set_' ts '.mat'], ...
    'mi',...
    'mcmc_info', 'data_info', 'mai', 'ri');  

marray = mcmc_get_walkers(tstamp,nIterID, projdir);
currmi = mi(2);
mvarray = masterVecArray(marray(:,:,1:100:end), mai);
marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,end);

parvals = struct('paramNames', currmi.namesOrd,...
    'paramVals', num2cell((marrayOrd(:, end, end))),...
    'reactionString', 'global');

% simulate model with that parameter point
Mobj = model_dsg2014_regen('plotmode', true, 'initialDNA', 30, 'paramInfo', parvals);

doseStruct = struct('dosedNames', {{'DNA p70--utr1--deGFP'}}, ...
    'dosedVals', 30);
%
% plot trajectories for that model . 
mcmc_plot_txtl(Mobj, parvals,doseStruct, ...
    'model_info', currmi, 'master_info', mai); 
