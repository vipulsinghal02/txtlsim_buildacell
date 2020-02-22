%% Prediction script for the IFFL.
% Vipul Singhal, Jun 2019
%
% % best way to do this is to set up a proper project, complete with dosing
% % etc. That will export the model for you too.
%
% % Generate model object
% mIFFL = model_txtl_lastetIFFL;
%
% % Get Data
% di = ZachIFFL_testdata;
%
%
%
% % set parameters that will be fixed. And the ones that will be sampled from
% % the training set E are labeled with "sample from training E".
%
% % the values tagged with verifiedAgainstTrainingE below are in agreement
% % with the values in trainingE and in mcmc_info_ZSIFFL_predictionA
%
% activeNames_used_to_construct_things = {...
% 'TX_elong_glob', exp(3.121)     % %                                         verifiedAgainstTrainingE
% 'TL_elong_glob', exp(3.436)     % %                                         verifiedAgainstTrainingE
% 'AGTPdeg_time', exp(10.05) % %                                              verifiedAgainstTrainingE
% 'AGTPreg_ON', exp( -3.9120) % %                                             verifiedAgainstTrainingE
% 'AGTPdeg_rate', exp(-9.7873) % % from the vnprl_F2 set,                     verifiedAgainstTrainingE
% 'TXTL_INDUCER_LASR_AHL_Kd', exp(13) % %                                     verifiedAgainstTrainingE
% 'TXTL_INDUCER_LASR_AHL_F', exp(0) % %                                       verifiedAgainstTrainingE
% 'TXTL_UTR_UTR1_Kd', exp(0.0542) % % from the vnprl_F2 set,                  verifiedAgainstTrainingE
% 'TXTL_PLAC_RNAPbound_Kd', exp(10.5) % sample from training E OR exp(10.5)   ??
% 'TXTL_PLAC_RNAPbound_F', exp( 1.5000) % %                                   verifiedAgainstTrainingE
% 'TXTL_NTP_RNAP_1_Kd', exp( 2.9459) % %                                      verifiedAgainstTrainingE
% 'TXTL_NTP_RNAP_2_Kd', exp( 13.9970) % %                                     verifiedAgainstTrainingE
% 'TL_AA_Kd', exp( 6.5566) % %                                                verifiedAgainstTrainingE
% 'TL_AGTP_Kd', exp( 14.5090) % %                                             verifiedAgainstTrainingE
% 'TXTL_RNAdeg_Kd', exp(15.6349) % % from the vnprl_F2 set                    verifiedAgainstTrainingE
% 'TXTL_INDUCER_TETR_ATC_Kd', exp(-2) % %                                     verifiedAgainstTrainingE
% 'TXTL_INDUCER_TETR_ATC_F', exp(1.577) % %                                   verifiedAgainstTrainingE
% 'TXTL_DIMER_tetR_Kd', exp(-10) % %                                          verifiedAgainstTrainingE
% 'TXTL_DIMER_tetR_F', exp(1.447) % %                                         verifiedAgainstTrainingE
% 'TXTL_PLAS_RNAPbound_Kd', exp(32) % sample from training E
% 'TXTL_PLAS_RNAPbound_F', exp(0) % %                                         verifiedAgainstTrainingE
% 'TXTL_PLAS_TFBIND_Kd', exp(6) % sample from training E
% 'TXTL_PLAS_TFRNAPbound_Kd', exp(9) % sample from training E
% 'TXTL_PROT_deGFP_MATURATION', exp( -6.0748) % %                             verifiedAgainstTrainingE
% 'TXTL_UTR_UTR1_F', exp( -0.2000) % %                                        verifiedAgainstTrainingE
% 'TXTL_PLAS_TFRNAPbound_F', exp(0) % %                                       verifiedAgainstTrainingE
% 'TXTL_PLAS_TFBIND_F', exp(0) % %                                            verifiedAgainstTrainingE
% 'TXTL_RNAPBOUND_TERMINATION_RATE', 0.15 % sample from training E
% 'TXTL_NTP_RNAP_1_F', exp(      0) % %                                       verifiedAgainstTrainingE
% 'TXTL_NTP_RNAP_2_F', exp(      0) % %                                       verifiedAgainstTrainingE
% 'TXTL_PTET_sequestration_Kd', exp(-0.5) % %                                 verifiedAgainstTrainingE
% 'TXTL_COMBINATORIAL_ACTIVATOR_KNOCKOFF_Kd', 0.1 % reverse rate here is actually the ASSOCIATION.                             % ??
% 'TXTL_PTET_sequestration_F', exp(1.314) % %                                 verifiedAgainstTrainingE
% 'TXTL_COMBINATORIAL_ACTIVATOR_KNOCKOFF_F', 1  % forward rate is the dissociation. want this to be 100 - 1000x higher than the                               % ??
% 'TL_AA_F', exp( -0.3000) % %                                                verifiedAgainstTrainingE
% 'TL_AGTP_F', exp( -1.2000) % %                                              verifiedAgainstTrainingE
% 'TXTL_RIBOBOUND_TERMINATION_RATE', 40 % sample from training E
% 'TXTL_RNAdeg_F', exp(0) % %                                                 verifiedAgainstTrainingE
% 'TXTL_RNAdeg_kc', exp(-0.2251) % % from the vnprl_F2 set                    verifiedAgainstTrainingE
% 'RNAP', exp(5.9) % sample from training E OR exp(5.9)                       verifiedAgainstTrainingE, but can sample
% 'RNase', exp(9.2) % %                                                       verifiedAgainstTrainingE
% 'Ribo', exp(5.9)... % sample from training E OR exp(9.2)                    verifiedAgainstTrainingE, but can sample
% };
%
% %% Export model and get it ready for simulation.
%
%
%
%
% %%
%
%     titls = arrayfun(@num2str, mi(1).dosedVals, 'UniformOutput', false);
%     mcmc_trajectories(mi(1).emo,...
%     di(mi(1).dataToMapTo),...
%     mi(1),...
%     (mod1(mi(1).paramMaps(mi(1).orderingIx,1))'),...
%     titls',...
%     {},...
%     'SimMode', 'curves')
% %
%
% titls = arrayfun(@num2str, mi(2).dosedVals, 'UniformOutput', false);
% titls_array = cell(length(titls), 1, length(mi(2).measuredSpeciesIndex));
% for i = 1:length(mi(2).measuredSpeciesIndex)
%     for j = 1:length(titls)
%         titls_array(j, 1, i) = titls(j);
%     end
% end
% mcmc_trajectories(mi(2).emo,...
%     di(mi(2).dataToMapTo),...
%     mi(2),...
%     (mod1(mi(2).orderingIx)'),...
%     titls_array,...
%     {},...
%     'SimMode', 'curves')
% % Set random seed, and sample model parameters randomly from the trainingE
% % dataset.
%
% % Generate trajectories from the sampled values, along with means and
% % standard deviations.
%
% % Generate trajectories in collated mode, without the standard deviations.
%
% %

%% Code construction:

% The best way to predict using the parameters estimated in the training E
% example is to set up an estimation problem with all the parameters is
% training E to be estimated, run it for a small amount of time, then
% instead of using the paramters estimated, use the parameters from
% training E.

% File and directory info:
%
% Project name:
%  'proj_ZSIFFL_predictionA'
%
% Directory where the project file is stored:
%  '/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/mcmc_simbio/projects'
%
% Directory where data will be stored:
%  '/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/mcmc_simbio/projects/proj_ZSIFFL_predictionA'
%
% Timestamp for this run (yyyymmdd_HHMMSS):
%  '20190613_131814'
%
% Project directory already exists, using this to store data
%  (in a subdirectory named 'simdata_20190613_131814_1_327').

% the training E project directory
close all
clear all
saveFinalFigs = '/Users/vipulsinghal/Dropbox/Documents/a_Journal_Papers/Drafts/txtl_bmc_bioinformatics/figs/Jul20_2019_v2/'
finafigmode = true
trainingEdir = [pwd '/mcmc_simbio/projects/proj_ZSIFFL_trainingE'];
projdir = [pwd '/mcmc_simbio/projects/proj_ZSIFFL_predictionA'];

% the array of parameters with all the fixed parameters and the
% parameters from the training E dataset.
addpath(projdir)
sls = regexp(projdir, '/', 'split');
extrastring = sls{end};
jpgsave = true;
figsave = false;

% Load model, mcmc_info, and data_info.
% construct simbiology model object(s)
mIFFL = model_txtl_lastetIFFL;
% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_predictionA(mIFFL);
di = ZachIFFL_testdata('all_trajectories');

tsIDtouse = 1;
plotflag = true;
switch tsIDtouse
    case 1 % this is after including the termination parameters.
        ts1 = '20190512_033129_1_1476';
        ts2 = '20190512_064547_1_1476';
        ts3 = '20190512_094712_1_1845';
        ts4 = '20190512_113712_1_1845';
        ts5 = '20190512_155207_1_1845';
        tstamp = {ts1 ts2 ts3 ts4 ts5};
        nIterID = {1:3 1:8 1:2 1:5 1:7};
        %         load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
        %             'mi',...
        %             'mcmc_info', 'data_info',  'ri');
end
tsToSave = '20190630_183715_1_327';
tsToSave = '20190630_183715_1_327';
tsToSave = '20190720_124011_1_327';
load(['/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/mcmc_simbio/projects/proj_ZSIFFL_predictionA'...
    '/simdata_' tsToSave '/full_variable_set_' tsToSave '.mat'], ...
    'mi',...
    'mcmc_info', 'data_info',  'ri');
data_info = di; % overwrite with new di, with all three trajectories, not just the mean. 
mai = mcmc_info.master_info;
% mi = mcmc_info.model_info
mai.masterVector
marray_full = mcmc_get_walkers(tstamp,nIterID, trainingEdir);
marray = marray_full(:,:,1:end);
clear marray_full
parnames = [...
    {'pol_{Kd, tet}'}           % 7
    {'pol_{Kd, lac}'}           % 21
    {'pol_{term}'}              % 23
    {'Ribo_{term}'}             % 28
    {'pol'}                     % 31
    {'Ribo'}                    % 33
    {'pol_{Kd,las}'}            % 37
    {'plas_{tf, Kd}'}           % 39
    {'plas-pol_{tf, Kd}'}    ]; % 40
% verified jun 13, 2019.

%%
%
% 'RNAP'
% 'Ribo'
% 'RecBCD'
% 'RNase'
% 'AGTP'
% 'CUTP'
% 'AA'
%
%
% 'protein lasR'
% 'OC12HSL:protein lasR'
%
%
% 'RNA utr1--lasR'
% 'Ribo:RNA utr1--lasR'
%
%
% 'DNA plac--utr1--lasR'
% 'RNAP:DNA plac--utr1--lasR'
% 'CUTP:AGTP:RNAP:DNA plac--utr1--lasR'
% 'term_RNAP:DNA plac--utr1--lasR'
%
%
% 'AA:AGTP:Ribo:RNA utr1--lasR'
% 'term_Ribo:RNA utr1--lasR'
%
%
% 'protein tetR'
% 'aTc'
% 'protein tetRdimer'
%
%
% 'RNA utr1--tetR'
% 'Ribo:RNA utr1--tetR'
%
%
% 'DNA plas--utr1--tetR'
% 'RNAP:DNA plas--utr1--tetR'
%
%
% 'AA:AGTP:Ribo:RNA utr1--tetR'
% 'term_Ribo:RNA utr1--tetR'
%
%
% 'protein deGFP'
% 'protein deGFP*'
%
%
% 'RNA utr1--deGFP'
% 'Ribo:RNA utr1--deGFP'
%
%
% 'DNA plas_ptet--utr1--deGFP'
% 'RNAP:DNA plas_ptet--utr1--deGFP'
%
%
% 'AA:AGTP:Ribo:RNA utr1--deGFP'
% 'term_Ribo:RNA utr1--deGFP'
%
%
% 'OC12HSL'
%
%
% 'AGTP:RNAP:DNA plac--utr1--lasR'
% 'CUTP:RNAP:DNA plac--utr1--lasR'
%
%
% 'AA:Ribo:RNA utr1--lasR'
%
%
% 'AGMP'
%
%
% 'RNA utr1--lasR:RNase'
% 'CUMP'
% 'Ribo:RNA utr1--lasR:RNase'
% 'AA:AGTP:Ribo:RNA utr1--lasR:RNase'
% 'term_Ribo:RNA utr1--lasR:RNase'
% 'AA:Ribo:RNA utr1--lasR:RNase'
%
%
% '2 aTc:protein tetRdimer'
%
%
% 'DNA plas--utr1--tetR:OC12HSL:protein lasR'
% 'RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR'
%
%
% 'CUTP:AGTP:RNAP:DNA plas--utr1--tetR'
% 'term_RNAP:DNA plas--utr1--tetR'
%
%
% 'CUTP:AGTP:RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR'
% 'term_RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR'
%
% 'AGTP:RNAP:DNA plas--utr1--tetR'
% 'CUTP:RNAP:DNA plas--utr1--tetR'
%
%
% 'AGTP:RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR'
% 'CUTP:RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR'
%
%
% 'AA:Ribo:RNA utr1--tetR'
%
%
% 'RNA utr1--tetR:RNase'
% 'Ribo:RNA utr1--tetR:RNase'
% 'AA:AGTP:Ribo:RNA utr1--tetR:RNase'
% 'term_Ribo:RNA utr1--tetR:RNase'
% 'AA:Ribo:RNA utr1--tetR:RNase'
%
%
% 'DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
% 'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
% 'CUTP:AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
% 'term_RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
% 'AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
% 'CUTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%
%
% 'DNA plas_ptet--utr1--deGFP:protein tetRdimer'
% 'DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer'
% 'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer'
%
%
% 'AA:Ribo:RNA utr1--deGFP'
%
%
% 'RNA utr1--deGFP:RNase'
% 'Ribo:RNA utr1--deGFP:RNase'
% 'AA:AGTP:Ribo:RNA utr1--deGFP:RNase'
% 'term_Ribo:RNA utr1--deGFP:RNase'
% 'AA:Ribo:RNA utr1--deGFP:RNase'

%%
mvarray = masterVecArray(marray(2:end,:,:), mai);
%     clear marray
for miID = 1:5 %1:length(mi)%1:
    currmi = mi(miID);
    dvStr = arrayfun(@num2str, currmi.dosedVals, 'UniformOutput', false);
    % titles array: # dose combs x (exp, sim) ==2 x # measured species.
    nms = length(currmi.measuredSpeciesIndex);
    ndc = size(currmi.dosedVals,2); %# dose combinations
    titls_array = cell(ndc, 2, nms);
    currdi = data_info(currmi.dataToMapTo);
    dn = currdi.dosedNames; %usually the last one will be the thing that is changing...
    % can generalize to use all dose info if it becomes needed.
    for msID = 1:nms
        % for each measured species (ms), plot the trajectories over
        % all the doses, for both experiment and simulation
        if iscell(data_info(currmi.dataToMapTo).measuredNames{msID})
            if ischar( data_info(currmi.dataToMapTo).measuredNames{msID}{1})
                ms = data_info(currmi.dataToMapTo).measuredNames{msID}{1};
            else
                error('what is your data type?')
            end
        elseif ischar( data_info(currmi.dataToMapTo).measuredNames{msID})
            ms = data_info(currmi.dataToMapTo).measuredNames{msID};
        else
            error('what is your data type?')
        end
        for dcID = 1:ndc
            dosestr = [];
            for dnID = 1:length(dn)
                dosestr = [dosestr ' ' dn{dnID} ' ' num2str(currdi.dosedVals(dnID, dcID))];
            end
            % dose combination by (exp, sim) by measured species ID
            titls_array{dcID, 1, msID} = ['Exp ' ms dosestr];
            titls_array{dcID, 2, msID} = ['Sim ' ms dosestr];
        end
    end
    titls_array
    samplePoints = ceil(size(mvarray, 3) * [.95, 1]);
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
    if finafigmode == false
        mcmc_trajectories(currmi.emo, currdi, currmi, marrayOrd,...
            titls_array, {},...
            'SimMode', 'mean', 'separateExpSim', true,'collateDoses', false,...
            'savematlabfig', figsave, 'savejpeg', jpgsave,...
            'projdir', projdir, 'tstamp', tsToSave,...
            'extrafignamestring', [num2str(miID) ' ' num2str(msID)]);%,'collateDoses', true,
        %     here we plot the individual trajectories for each species.
       
    end
     ms = currmi.measuredSpecies;
    
end


%%
% 
% ms = {{'protein lasR'};
%     {'OC12HSL:protein lasR'};
%     {'OC12HSL'};
%     {'RNA utr1--lasR';
%     'Ribo:RNA utr1--lasR';
%     'AA:Ribo:RNA utr1--lasR';
%     'AA:AGTP:Ribo:RNA utr1--lasR';
%     'term_Ribo:RNA utr1--lasR'};
%     {'protein tetR'};
%     {'protein tetRdimer'};
%     {'DNA plas_ptet--utr1--deGFP:protein tetRdimer';
%     'DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer';
%     'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer'};
%     {'DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'CUTP:AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'term_RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'CUTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'};
%     {'protein deGFP*', 'protein deGFP'}};
% 
% for miID = 1:5
%     currmi = mi(miID);
%     currdi = data_info(currmi.dataToMapTo);
%     tv = currdi.timeVector;
%     dose  = currmi.dosedVals';
%     [da{miID}, idxnotused{miID}] = simulatecurves(currmi.emo,marrayOrd(:,:)', 20, dose, tv, ms);
% end
% 
% % {'2 aTc:protein tetRdimer'};
% % {'aTc'};
% 
% close all
% currda = da{4};
% specieslabels = {'lasR', 'AHL:lasR','AHL', 'lasR RNA', 'tetR', 'tetRdimer',...
%     'repressed_comb_prom', 'activated_comb_prom', 'GFP'}; % 'aTc:tetR','aTc',
% doselabels = {'1nM', '0.1nM', '0.01nM', '0.001nM', '0.0001nM', '0.00001nM', '0nM'};
% figure
% nSp = size(currda, 2);
% [fac1, fac2] = twofactors(nSp);
% 
% for speciesGrp = 1:nSp
%     subplot(fac1,  fac2, speciesGrp)
%     for dID = 1:length(doselabels)
%         plot(tv, mean(currda(:, speciesGrp, :, dID), 3) + 0.01*randn(length(tv),1), 'LineWidth', 1.5)
%         hold on
%     end
%     title([specieslabels{speciesGrp}])
%     legend(doselabels)
% end
% 
% %% same thing, explore the lasR activation dynamics.
% % why is there significant activation even at 0.03nM lasR.
% %should it not just be close to 0 expression?
% 
% 
% ms = {{'protein lasR'};
%     {'OC12HSL:protein lasR'};
%     {'OC12HSL'};
%     {'RNA utr1--lasR';
%     'Ribo:RNA utr1--lasR';
%     'AA:Ribo:RNA utr1--lasR';
%     'AA:AGTP:Ribo:RNA utr1--lasR';
%     'term_Ribo:RNA utr1--lasR'};
%     {'protein tetR'};
%     {'protein tetRdimer'};
%     {'DNA plas_ptet--utr1--deGFP:protein tetRdimer';
%     'DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer';
%     'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer'};
%     {'DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'CUTP:AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'term_RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'CUTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'};
%     {'protein deGFP*', 'protein deGFP'}};
% 
% % {'2 aTc:protein tetRdimer'};
% % {'aTc'};
% for miID = 1:5
%     
%     currmi = mi(miID);
%     currdi = data_info(currmi.dataToMapTo);
%     tv = currdi.timeVector;
%     dose  = currmi.dosedVals';
%     if miID == 2
%         dose(end, 1) = 0
%     end
%     
%     [da{miID}, idxnotused{miID}] = simulatecurves(currmi.emo,marrayOrd(:,:)', 20, dose, tv, ms);
% end
% close all
%%
% currda = da{4};
% specieslabels = {'lasR', 'AHL:lasR','AHL', 'lasR RNA', 'tetR', 'tetRdimer',...
%     'repressed_comb_prom', 'activated_comb_prom', 'GFP'}; % 'aTc:tetR','aTc',
% doselabels = {'2nM', '1nM', '0.5nM', '0.25nM', '0.125nM', '0.0625nM', '0.03125nM'};
% figure
% nSp = size(currda, 2);
% [fac1, fac2] = twofactors(nSp);
% 
% for speciesGrp = 1:nSp
%     subplot(fac1,  fac2, speciesGrp)
%     for dID = 1:length(doselabels)
%         plot(tv, mean(currda(:, speciesGrp, :, dID), 3) + 0.01*randn(length(tv),1), 'LineWidth', 1.5)
%         hold on
%     end
%     title([specieslabels{speciesGrp}])
%     legend(doselabels)
% end

% I dont think this can be fixed. Maybe if we had the lasR DNA perturbation
% data. Don't have that. Not gonna sweat it.

%% Generate prediction figure
% Next, We plot the collated results: two columns, one for data, the other
% for the predictions.

% resimulate, but with just the matured degfp.
ms = {{'protein deGFP*'}};
NN = 50; % number of trajectories.
da = cell(NN, 1);
idxnotused = cell(NN, 1);
for miID = 1:5
    currmi = mi(miID);
    currdi = data_info(currmi.dataToMapTo);
    tv = currdi.timeVector;
    dose  = currmi.dosedVals';
    if miID == 2
        dose(end, 1) = 0;
    end
    [da{miID}, idxnotused{miID}] = simulatecurves(currmi.emo,marrayOrd(:,:)', NN, dose, tv, ms);
end
%% plot the trajectories figure. 
close all


%
ylims = {[4, 4, 4, 4, 12];
    2*[4, 4, 4, 4, 12];
    3*[4, 4, 4, 4, 12]};


lengthToPlotArray = [31, 41, 81];
for outercount = 1:length(lengthToPlotArray)
%     outercount = outercount+1;
    lengthToPlot = lengthToPlotArray(outercount)
        
    figure
    ss = get(0, 'screensize');
    set(gcf, 'Position', [ss(3)*(1-1/1.3) ss(4)*(1-1/1.3) ss(3)/3.5 ss(4)/1.1]);
    
    miToUse = [1 2 3 4 5];
    titleArray = {'3OC12';
        'lasR DNA';
        'aTc';
        'tetR DNA';
        'deGFP DNA'};
    
    legends = [{'10uM', '1uM', '.1uM', '10nM', '1nM', '.1nM', '0nM'};
        {'2nM', '1nM', '.5nM', '.25nM', '125pM', '62.5pM', '0pM'};
        {'10uM', '1uM', '.1uM', '10nM', '1nM', '0.1nM', '0nM'};
        {'1nM', '.1nM', '.01nM', '1pM', '.1pM', '10fM', '0fM'};
        {'4nM', '2nM', '1nM', '.5nM', '.25nM', '.125nM', '0nM'}];
    ms = currmi.measuredSpecies;
    mvarray = masterVecArray(marray(2:end,:,:), mai);
%     lengthToPlot = 31; % 180 min = 31
    %     clear marray
    colorz = flipud(parula(length(legends(1, :))+2));
    for count = 1:length(miToUse)  %1:length(mi)%1:
        currmi = mi(miToUse(count));
        currdi = data_info(currmi.dataToMapTo);
        %     samplePoints = ceil(size(mvarray, 3) * [.95, 1]);
        %     marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
        tv = currdi.timeVector;
        dose  = currmi.dosedVals';
        if miToUse(count) == 2
            dose(end, 1) = 0;
        end
        %[da{count}, idxnotused{count}] = simulatecurves(currmi.emo,marrayOrd(:,:)', 30, dose, tv, ms);
        
        subplot(5, 2, (count-1)*2+1)
        for dID = 1:length(legends(count, :))
            plot(tv(1:lengthToPlot)/60, 0.001*mean(currdi.dataArray(1:lengthToPlot, 1, :, dID),3), 'LineWidth', 1.5, 'Color', colorz(dID+2, :))
            hold on
        end
        title([titleArray{count} ' varying (Experiment)'])
        ylabel('deGFP, uM', 'FontSize', 14)
        axis([0 (lengthToPlot-1)*6, 0 ylims{outercount}(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
        
        subplot(5, 2, (count-1)*2+2)
        for dID = 1:length(legends(count, :))
            dID;
            plot(tv(1:lengthToPlot)/60, 0.001*mean(da{count}(1:lengthToPlot, 1, :, dID), 3), 'LineWidth', 1.5, 'Color', colorz(dID+2, :))
            hold on
        end
        title([titleArray{count} ' varying (Prediction)'])
        legend(legends(count, :), 'Location', 'NorthWest', 'FontSize', 13)
        legend('boxoff')
        
        %     ylabel('deGFP, nM')
        axis([0 (lengthToPlot-1)*6 0 ylims{outercount}(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
    end
    
%     print([saveFinalFigs 'prediction_traj_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-depsc')
%     print([saveFinalFigs 'prediction_traj_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-djpeg')
%     
%     print([saveFinalFigs 'prediction_traj_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-dpng')
    %
    
    %
    % Now generate the plots
    %
    % build the data array
    
    figure
    % create the endpoint plots
    % The first location is empty, for the IFFL schematic.
    exp_endpoints = zeros(size(legends)); % number of experiments by number of doses.
    sim_endpoints = zeros(size(legends));
    sim_endpointsSD = zeros(size(legends));
    exp_endpointsSD = zeros(size(legends));
    doseArray = zeros(size(legends));
    for i = 1:5
        doseArray(i, :) = data_info(i).dosedVals;
    end
    doseArray(1, end) = 0.1/100; % 30c12
    doseArray(2, end) = 0.0625/4; % lasR dna
    doseArray(3, end) = 0.1/100; % aTc
    doseArray(4, end) = 0.00001/100; % tetR dna
    doseArray(5, end) = 0.125/4; % deGFP dna
    
    
    for count = 1:length(miToUse)
        currmi = mi(miToUse(count));
        currdi = data_info(currmi.dataToMapTo);
        
        % get the experimental data endpoints
        for dID = 1:length(legends(count, :))
            exp_endpoints(count, dID) = ...
                0.001*mean(currdi.dataArray(lengthToPlot, 1, :, dID),3);
            exp_endpointsSD(count, dID) = ...
               std(0.001*currdi.dataArray(lengthToPlot, 1, :, dID),0,3);
            sim_endpoints(count, dID) = ...
                0.001*mean(da{count}(lengthToPlot, 1, :, dID), 3, 'omitnan');
            sim_endpointsSD(count, dID) = ...
                std(0.001*da{count}(lengthToPlot, 1, :, dID),0,  3, 'omitnan');
        end
    end
    
    
    
    figure
    ss = get(0, 'screensize');
    set(gcf, 'Position', [ss(3)*(1-1/1.3) ss(4)*(1-1/1.3) ss(3)/2.5 ss(4)/2]);
    
    
    xTickLabels = [{'10', '10^{-1}', '10^{-3}', '0'};
        {'2', '0.5', '0.125', '0'};
        {'10',  '10^{-1}',  '10^{-3}',  '0'};
        {'1', '10^{-2}',  '10^{-4}',  '0'};
        {'4', '1', '2^{-2}',  '0'}];
    titleArray = {'3OC12';
        'lasR DNA';
        'aTc';
        'tetR DNA';
        'deGFP DNA'};
    xLab = {'3OC12, uM';
        'lasR DNA, nM';
        'aTc, uM';
        'tetR DNA, nM';
        'deGFP DNA, nM'};
    
    for count = 1:length(miToUse)
        
        subplot(2, 3, (count+1))
        ax = gca
        ax.XScale = 'log';
        errorbar(doseArray(count, :),...
            exp_endpoints(count, :),...
            exp_endpointsSD(count, :),...
            'LineWidth', 2);
        
        hold on
        
        errorbar(doseArray(count, :), ...
            sim_endpoints(count, :), ...
            sim_endpointsSD(count, :), ...
            'LineWidth', 2)
        ax = gca
        ax.XScale = 'log';
        
        
        
%         
%         ax = semilogx(doseArray(count, :), ...
%             exp_endpoints(count, :)/1000, ...
%             'LineWidth', 2);
%         hold on
%         errorbar(doseArray(count, :), ...
%             sim_endpoints(count, :), ...
%             sim_endpointsSD(count, :), 'LineWidth', 2)
        hold on
        if count ==1
            legend('Experiment', 'Prediction', 'Location', 'NorthWest', ...
                'FontSize', 18)
        end
        
        ax = gca
        ax.XTick = fliplr(doseArray(count, 1:2:end));
        ax.XTickLabel = fliplr(xTickLabels(count, :));
        ax.FontSize = 18
        %     title(titleArray{count})
        ylabel('deGFP, uM');
        xlabel(xLab{count});
        
    end
    print([saveFinalFigs 'prediction_end_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-depsc')
%     print([saveFinalFigs 'prediction_end_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-djpeg')
%     print([saveFinalFigs 'prediction_end_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-dpng')
end

close all

%%

%
%
%
%
%     % get the simulated data endpoints.
%     title([titleArray{count} ' varying (Experiment)'])
%     ylabel('deGFP, uM', 'FontSize', 14)
%     axis([0 (lengthToPlot-1)*6, 0 ylims(count)])
%     if count ==length(miToUse)
%         xlabel('time, minutes', 'FontSize', 14)
%     end
%     ax = gca;
%     ax.FontSize = 14;
%
%     subplot(5, 2, (count-1)*2+2)
%     for dID = 1:length(legends(count, :))
%         dID;
%         plot(tv(1:lengthToPlot)/60, 0.001*mean(da{count}(1:lengthToPlot, 1, :, dID), 3), 'LineWidth', 1.5, 'Color', colorz(dID+2, :))
%         hold on
%     end
%     title([titleArray{count} ' varying (Prediction)'])
%     legend(legends(count, :), 'Location', 'NorthWest', 'FontSize', 13)
%     legend('boxoff')
%
%     %     ylabel('deGFP, nM')
%     axis([0 (lengthToPlot-1)*6 0 ylims(count)])
%     if count ==length(miToUse)
%         xlabel('time, minutes', 'FontSize', 14)
%     end
%     ax = gca;
%     ax.FontSize = 14;
%
%
%
%




%% Simulate and plot all the species in the system that could be of interest.
% The geometry that is failing is, for example, number 3, 4.
%    SimBiology Species Array
%
%    Index:    Compartment:    Name:                                                                     InitialAmount:    InitialAmountUnits:
%    1         contents        RNAP                                                                      100
%    2         contents        Ribo                                                                      30
%    3         contents        RecBCD                                                                    5
%    4         contents        RNase                                                                     100
%    5         contents        AGTP                                                                      3.18005e+06
%    6         contents        CUTP                                                                      1.90803e+06
%    7         contents        AA                                                                        3.18005e+07
%    8         contents        protein lasR                                                              0
%    9         contents        OC12HSL:protein lasR                                                      0
%    10        contents        RNA utr1--lasR                                                            0
%    11        contents        Ribo:RNA utr1--lasR                                                       0
%    12        contents        DNA plac--utr1--lasR                                                      0
%    13        contents        RNAP:DNA plac--utr1--lasR                                                 0
%    14        contents        CUTP:AGTP:RNAP:DNA plac--utr1--lasR                                       0
%    15        contents        term_RNAP:DNA plac--utr1--lasR                                            0
%    16        contents        AA:AGTP:Ribo:RNA utr1--lasR                                               0
%    17        contents        term_Ribo:RNA utr1--lasR                                                  0
%    18        contents        protein tetR                                                              0
%    19        contents        aTc                                                                       0
%    20        contents        protein tetRdimer                                                         0
%    21        contents        RNA utr1--tetR                                                            0
%    22        contents        Ribo:RNA utr1--tetR                                                       0
%    23        contents        DNA plas--utr1--tetR                                                      0
%    24        contents        RNAP:DNA plas--utr1--tetR                                                 0
%    25        contents        AA:AGTP:Ribo:RNA utr1--tetR                                               0
%    26        contents        term_Ribo:RNA utr1--tetR                                                  0
%    27        contents        protein deGFP                                                             0
%    28        contents        protein deGFP*                                                            0
%    29        contents        RNA utr1--deGFP                                                           0
%    30        contents        Ribo:RNA utr1--deGFP                                                      0
%    31        contents        DNA plas_ptet--utr1--deGFP                                                0
%    32        contents        RNAP:DNA plas_ptet--utr1--deGFP                                           0
%    33        contents        AA:AGTP:Ribo:RNA utr1--deGFP                                              0
%    34        contents        term_Ribo:RNA utr1--deGFP                                                 0
%    35        contents        OC12HSL                                                                   0
%    36        contents        AGTP:RNAP:DNA plac--utr1--lasR                                            0
%    37        contents        CUTP:RNAP:DNA plac--utr1--lasR                                            0
%    38        contents        AA:Ribo:RNA utr1--lasR                                                    0
%    39        contents        AGMP                                                                      0
%    40        contents        RNA utr1--lasR:RNase                                                      0
%    41        contents        CUMP                                                                      0
%    42        contents        Ribo:RNA utr1--lasR:RNase                                                 0
%    43        contents        AA:AGTP:Ribo:RNA utr1--lasR:RNase                                         0
%    44        contents        term_Ribo:RNA utr1--lasR:RNase                                            0
%    45        contents        AA:Ribo:RNA utr1--lasR:RNase                                              0
%    46        contents        2 aTc:protein tetRdimer                                                   0
%    47        contents        DNA plas--utr1--tetR:OC12HSL:protein lasR                                 0
%    48        contents        RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR                            0
%    49        contents        CUTP:AGTP:RNAP:DNA plas--utr1--tetR                                       0
%    50        contents        term_RNAP:DNA plas--utr1--tetR                                            0
%    51        contents        CUTP:AGTP:RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR                  0
%    52        contents        term_RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR                       0
%    53        contents        AGTP:RNAP:DNA plas--utr1--tetR                                            0
%    54        contents        CUTP:RNAP:DNA plas--utr1--tetR                                            0
%    55        contents        AGTP:RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR                       0
%    56        contents        CUTP:RNAP:DNA plas--utr1--tetR:OC12HSL:protein lasR                       0
%    57        contents        AA:Ribo:RNA utr1--tetR                                                    0
%    58        contents        RNA utr1--tetR:RNase                                                      0
%    59        contents        Ribo:RNA utr1--tetR:RNase                                                 0
%    60        contents        AA:AGTP:Ribo:RNA utr1--tetR:RNase                                         0
%    61        contents        term_Ribo:RNA utr1--tetR:RNase                                            0
%    62        contents        AA:Ribo:RNA utr1--tetR:RNase                                              0
%    63        contents        DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR                           0
%    64        contents        RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR                      0
%    65        contents        CUTP:AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR            0
%    66        contents        term_RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR                 0
%    67        contents        AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR                 0
%    68        contents        CUTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR                 0
%    69        contents        DNA plas_ptet--utr1--deGFP:protein tetRdimer                              0
%    70        contents        DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer         0
%    71        contents        RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer    0
%    72        contents        AA:Ribo:RNA utr1--deGFP                                                   0
%    73        contents        RNA utr1--deGFP:RNase                                                     0
%    74        contents        Ribo:RNA utr1--deGFP:RNase                                                0
%    75        contents        AA:AGTP:Ribo:RNA utr1--deGFP:RNase                                        0
%    76        contents        term_Ribo:RNA utr1--deGFP:RNase                                           0
%    77        contents        AA:Ribo:RNA utr1--deGFP:RNase


