% Generate Incoherent Feed-Forward Loop characterization plots from the paper:
% 
% A MATLAB Toolbox for Modeling Genetic Circuits in Cell-Free Systems
% by Vipul Singhal, Zoltan A. Tuza, Zachary S. Sun, and Richard M. Murray
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. check that the working directory is the one where folders like core,
% components, mcmc_simbio, config etc are. This is the trunk directory. 
% 2. Download dataset, and place in the correct directory. You need around 800mb
% of free hard disk space to download the full dataset. Download at: 
% https://www.dropbox.com/sh/9zn7666ll6fm4po/AAB1nS2Vi3Aej7vmPD-kQBhXa?dl=0
% 3. Make sure you have 7 GB of free ram. Otherwise, set the flag 
% 'small_ram' to true, and you will plot the same things but with only
% 4 GB of ram needed. This only plots the last few iterations of the MCMC.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
small_ram = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set directories and load data and files. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txtl_init
mcmc_init

saveFinalFigs = [pwd '\generate_papaer_figs\charac']
projdir = [pwd '\mcmc_simbio\projects\proj_ZSIFFL_trainingC_v6'];
% Set working directory to the txtlsim toolbox directory.
addpath(projdir)
sls = regexp(projdir, '\', 'split');
extrastring = sls{end};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load model, mcmc_info, and data_info (experimental data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct simbiology model object(s)
mtet = model_txtl_ptetdeGFP_pLactetR_aTc;
mlas = model_txtl_pLacLasR_pLasdeGFP;
mlac = model_txtl_pLacdeGFP;
% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_training_fullC_v6(mtet, mlac, mlas);
di = data_ZSIFFL;
mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

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
            tstamp = {ts47 ts48 ts49 ts50 ts51 ts52 ts53 ts54 ts55 ts56 ...
                 ts57 ts58 ts59};
            nIterID = {1:40 1:29 1:40 1:40 1:40 1:40 1:100 1:21 1:100 1:60 ...
                1:100 1:100 1:3};  
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
        load([projdir '\simdata_' tstamp{end} '\full_variable_set_'...
            tstamp{end} '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info',  'ri');
end
tsToSave = tstamp{end};
mai.masterVector
marray = mcmc_get_walkers(tstamp,nIterID, projdir);

parnames = [...
    {'TX_{cat}'}
    {'TL_{cat}'}
    {'\tau_{atp}'}
    {'pol_{Kd, tet}'}
    {'rep_{Kd, tet}'}
    {'aTc_{Kd}'}
    {'pol_{Kd, lac}'}
    {'pol_{term}'}
    {'Ribo_{term}'}
    {'pol'}
    {'Ribo'}
    {'3OC12_{Kd}'}
    {'pol_{Kd,las}'}
    {'plas_{tf, Kd}'}
    {'plas-pol_{tf, Kd}'}    ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the trace plot and the posterior (marginal) plot
mcmc_plot(marray(:, 1:end,1:end), parnames(:),...
    'savematlabfig', false,'plotDistribution', false,...
    'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% marray is a # parameters x # walkers x # mcmc samples array. You can
% subset this for plotting, especially when plotting the posterior distribution 
%(which needs a lot of ram), like so: 
mcmc_plot(marray(:, 1:end,1:100:end), parnames(:),...
    'plotDistribution', true,...
    'savematlabfig', false, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mvarray = masterVecArray(marray, mai);
for miID = 1:length(mi)%1:
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
                dosestr = ...
                    [dosestr ' ' dn{dnID} ' '...
                    num2str(currdi.dosedVals(dnID, dcID))];
            end
            % dose combination by (exp, sim) by measured species ID
            titls_array{dcID, 1, msID} = ['Exp ' ms dosestr];
            titls_array{dcID, 2, msID} = ['Sim ' ms dosestr];
        end
    end
    samplePoints = ceil(size(mvarray, 3) * [1]);
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the endpoint and trajectory figs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ms = {{'protein deGFP*'}};
for miID = 1:5
    currmi = mi(miID);
    currdi = data_info(currmi.dataToMapTo);
    tv = currdi.timeVector;
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
    dose  = currmi.dosedVals';
    vi = get(currmi.emo, 'ValueInfo');
    for i = 1:length(vi)
        vi(i)
    end
    [da{miID}, idxnotused{miID}] = ...
        simulatecurves(currmi.emo,marrayOrd(:,:)', 50, dose, tv, ms);
end
ylims = {[15 5 5 10 15]};
lengthToPlotArray = [61];
for outercount = 1:length(lengthToPlotArray)
    lengthToPlot = lengthToPlotArray(outercount);
    set(0,'Units','normalized')
    figure
    set(gcf,'Units', 'normalized')
    set(gcf, 'Position', [0.05, 0.05, 0.2, 0.9])
    miToUse = [1 2 3 4 5];
    
    titleArray = {'pTet-deGFP DNA';
        'pLac-TetR DNA';
        'aTc';
        'pLac-deGFP DNA';
        '3OC12HSL'};
    
    legends = [{'4nM', '2nM', '1nM', '0.5nM', '250pM', '125pM', '0nM'}
        {'2nM', '.2nM', '.02nM', '2pM', '0.2pM', '20fM', '0nM'};
        {'10uM', '1uM', '0.1uM', '10nM', '1nM', '0.1nM', '0nM'}
        {'2nM', '1nM', '.5nM', '.25nM', '.125nM', '62.5pM', '0nM'};
        {'10uM', '1uM', '0.1uM', '10nM', '1nM', '0.1nM', '0nM'}];
    colorz = flipud(parula(length(legends(1, :))+2));
    for count = 1:length(miToUse)  %1:length(mi)%1:
        currmi = mi(miToUse(count));
        currdi = data_info(currmi.dataToMapTo);
        tv = currdi.timeVector;
        subplot(5, 2, (count-1)*2+1)
        for dID = 1:length(legends(count, :))
            plot(tv(1:lengthToPlot)/60,...
                0.001*mean(currdi.dataArray(1:lengthToPlot, 1, 1, dID),3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            
            hold on
        end
        title([titleArray{count} ' (Exp)'])
        ylabel('deGFP, uM', 'FontSize', 14)
        axis([0 (lengthToPlot-1)*8, 0 ylims{outercount}(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
        subplot(5, 2, (count-1)*2+2)
        for dID = 1:length(legends(count, :))
            
            plot(tv(1:lengthToPlot)/60,...
                0.001*mean(da{count}(1:lengthToPlot, 1, :, dID), 3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            hold on
        end
        title([titleArray{count} ' (Fit)'])
        legend(legends(count, :), 'Location', 'NorthWest', 'FontSize', 10)
        legend('boxoff')
        
        axis([0 (lengthToPlot-1)*8 0 ylims{outercount}(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
    end
    % create the endpoint plots
    % The first location is empty, for the IFFL schematic.
    figure
    exp_endpoints = zeros(size(legends)); % number of experiments by number of doses.
    exp_endpointsSD = zeros(size(legends));
    sim_endpoints = zeros(size(legends));
    sim_endpointsSD = zeros(size(legends));
    doseArray = zeros(size(legends));
    for i = 1:5
        doseArray(i, :) = data_info(i).dosedVals;
    end
    % doseArray(1, end) = 0.1/100; % 30c12
    % doseArray(2, end) = 0.0625/4; % lasR dna
    % doseArray(3, end) = 0.1/100; % aTc
    % doseArray(4, end) = 0.00001/100; % tetR dna
    % doseArray(5, end) = 0.125/4; % deGFP dna
    % build the data array
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
    
    xTickLabels = [{'4', '1', '2^{-2}',  '2^{-4}'};
        {'2', '.02',  '.0002',  '0'};
        {'10', '10^{-1}', '10^{-3}', '0'};
        {'2', '0.5', '0.125', '0'};
        {'10',  '10^{-1}',  '10^{-3}',  '0'}];
    titleArray = {'pTet-deGFP DNA (nM)';
        'pLac-tetR DNA (nM)';
        'aTc (uM)';
        'pLac-deGFP DNA (nM)';
        '3OC12 (uM)'};
    xLab = {'pTet-deGFP DNA (nM)';
        'pLac-tetR DNA (nM)';
        'aTc (uM)';
        'pLac-deGFP DNA (nM)';
        '3OC12 (uM)'};
    xLims = [0.0625            4
        2e-06 2
        0.01 10000
        0.03125 2
        0.01 10000];
    doseArray =...
        [   4            2            1          0.5         0.25        0.125       0.0625
            2          0.2         0.02        0.002       0.0002        2e-05            2e-06
        10000         1000          100           10            1          0.1            0.01
            2            1          0.5         0.25        0.125       0.0625      0.03125
        10000         1000          100           10            1          0.1            0.01];
    for count = 1:length(miToUse)
        subplot(2, 3, (count+1))
        ax = gca;
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
        ax = gca;
        ax.XScale = 'log';
        ax.XLim = xLims(count, :);
        hold on
        ax = gca;
        ax.XTick = fliplr(doseArray(count, 1:2:end));
        ax.XTickLabel = fliplr(xTickLabels(count, :));
        ax.FontSize = 16;
        %     title(titleArray{count})
        ylabel('deGFP, uM');
        xlabel(xLab{count});
    end
end
