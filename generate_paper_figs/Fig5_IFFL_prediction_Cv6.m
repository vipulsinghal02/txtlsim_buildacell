% Generate Incoherent Feed-Forward Loop prediction plots from the paper:
% 
% A MATLAB Toolbox for Modeling Genetic Circuits in Cell-Free Systems
% by Vipul Singhal, Zoltan A. Tuza, Zachary S. Sun, and Richard M. Murray
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Check that the working directory is the one where folders like core,
% components, mcmc_simbio, config etc are. This is the trunk directory. 
% 2. Download datasets, and place in the mcmc_simbio/projects directory. 
% You need around 800mb
% of free hard disk space to download the full dataset. Download at: 
% https://www.dropbox.com/sh/9zn7666ll6fm4po/AAB1nS2Vi3Aej7vmPD-kQBhXa?dl=0
% and 
% https://www.dropbox.com/sh/sycgp5dr2arxop0/AACjd9owdihvQEQpQwdHFb-ya?dl=0
% 3. Make sure you have 3 GB of free ram. Otherwise, set the flag 
% 'small_ram' to true, and you will plot the same things but with only
% 0.8GB of ram needed. This only plots the last few iterations of the MCMC.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set flags
small_ram = false;
txtl_init
mcmc_init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These two datasets must be downloaded from the dropbox links above, and 
% placed in the mcmc_simbio/projects directory
trainingEdir = [pwd '\mcmc_simbio\projects\proj_ZSIFFL_trainingC_v6'];
projdir = [pwd '\mcmc_simbio\projects\proj_ZSIFFL_predictionA'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the array of parameters with all the fixed parameters and the
% parameters from the training E dataset.
addpath(projdir)
sls = regexp(projdir, '\', 'split');
extrastring = sls{end};
% Load model, mcmc_info, and data_info.
% construct simbiology model object(s)
mIFFL = model_txtl_lastetIFFL;
% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_predictionA(mIFFL);
di = ZachIFFL_testdata('all_trajectories');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load simulation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsIDtouse = 1;
switch tsIDtouse
    case 1 
        ts31 = '20200612_214710_1_1476';
        ts32 = '20200614_005043_1_1476';
        ts33 = '20200614_005043_2_2214';
        ts34 = '20200614_005043_3_2953';
        ts35 = '20200614_005043_4_3691';
        ts36 = '20200614_005043_5_2953';
        ts37 = '20200615_233435_1_2953';
        ts38 = '20200616_183409_1_2953';
        ts39 = '20200616_183409_2_2214';
        ts40 = '20200619_213706_1_2214';
        ts41 = '20200619_213706_2_1476';
        ts42 = '20200621_100025_1_4429';
        ts43 = '20200622_121135_1_4429';
        ts44 = '20200623_175503_1_4429';
        ts45 = '20200627_132326_1_3691';
        ts46 = '20200627_132326_2_2953';
        ts47 = '20200627_132326_3_2583';
        ts48 = '20200627_132326_4_2214';
        ts49 = '20200703_000710_1_2030';
        ts50 = '20200703_000710_2_1845';
        ts51 = '20200703_000710_3_1661';
        ts52 = '20200703_000710_4_1476';
        ts53 = '20200710_230144_1_1476';
        ts54 = '20200710_230144_2_1292';
        ts55 = '20200716_111111_1_1292';
        ts56 = '20200716_111111_2_1107';
        ts57 = '20200723_011054_1_1107';
        ts58 = '20200723_011054_2_923';
        ts59 = '20200723_011054_3_738';
        if small_ram
            tstamp = {ts57 ts58 ts59};
            nIterID = {1:100 1:100 1:3};  
        else
            tstamp = ...
                {ts31 ts32 ts33 ts34 ts35 ts36...
                 ts37 ts38 ts39 ts40 ts41 ts42 ts43 ts44 ts45 ts46 ...
                 ts47 ts48 ts49 ts50 ts51 ts52 ts53 ts54 ts55 ts56 ...
                 ts57 ts58 ts59};
            nIterID = ...
                {1:3 1:10 1:10 1:10 1:10 1:5 ...
                 1:22 1:40 1:39 1:40 1:2 1:10 1:36 1:100 1:40 1:40 ...
                 1:40 1:29 1:40 1:40 1:40 1:40 1:100 1:21 1:100 1:60 ...
                 1:100 1:100 1:68};  
        end
end


tsToSave = '20200525_011005_1_327'; % DO NOT CHANGE. this is where the prediction data is kept.
load([pwd '\mcmc_simbio\projects\proj_ZSIFFL_predictionCv3'...
    '\simdata_' tsToSave '\full_variable_set_' tsToSave '.mat'], ...
    'mi',...
    'mcmc_info',  'ri');
data_info = di; % overwrite with new di, with all three trajectories, not just the mean. 
mai = mcmc_info.master_info;
marray = mcmc_get_walkers(tstamp,nIterID, trainingEdir);

% 
%%
mvarray = masterVecArray(marray([1:3 5:15], :,:), mai);
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
    samplePoints = ceil(size(mvarray, 3) * [.95, 1]);
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
     ms = currmi.measuredSpecies;
end
%%
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
%%

%
ylims = {3*[4, 4, 4, 4, 12]};
lengthToPlotArray = [81];
for outercount = 1:length(lengthToPlotArray)
%     outercount = outercount+1;
    lengthToPlot = lengthToPlotArray(outercount);
    set(0,'Units','normalized')
    figure
    set(gcf,'Units', 'normalized')
    set(gcf, 'Position', [0.05, 0.05, 0.2, 0.9])
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
    mvarray = masterVecArray(marray([1:3 5:15],:,:), mai);
    colorz = flipud(parula(length(legends(1, :))+2));
    for count = 1:length(miToUse)  %1:length(mi)%1:
        currmi = mi(miToUse(count));
        currdi = data_info(currmi.dataToMapTo);
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
    % UNCOMMEENT TO SAVE
	% print([saveFinalFigs 'prediction_traj_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-djpeg')

    figure
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

    set(0,'Units','normalized')
    figure
    set(gcf,'Units', 'normalized')
    set(gcf, 'Position', [0.05, 0.1, 0.35, 0.5])
    
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
    % UNCOMMEENT TO SAVE
    %print([saveFinalFigs 'prediction_endpoint_N' num2str(NN) '_' num2str((lengthToPlot-1)*6) 'min'],'-djpeg')
end
