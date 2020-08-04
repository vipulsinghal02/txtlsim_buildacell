

clear all
close all
saveFinalFigs = 'D:\Dropbox\Documents\a_Journal_Papers\Drafts\txtl_bmc_bioinformatics\figs\trainingC_v6\'
finafigmode = true
% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '\mcmc_simbio\projects\proj_ZSIFFL_trainingC_v6'];
addpath(projdir)
sls = regexp(projdir, '\', 'split');
extrastring = sls{end};
jpgsave = false;
figsave = false;

% Load model, mcmc_info, and data_info.
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
% plot data from existing simulations.


% 
% ts1 = '20200518_175057_1_7381';
% nIterID1 = {1:3};
% marray1 = mcmc_get_walkers({ts1},nIterID1, projdir);
% marray1 = marray1(:,1:1000, :);

tsIDtouse = 1;
plotflag = true;

switch tsIDtouse
    case 1 % this is after including the termination parameters.
        
        ts2 = '20200521_154748_1_73814';%100, 100
        ts3 = '20200521_154748_2_36907'; %100, 200 % mostly random walking until here
        ts4 = '20200521_154748_3_14763'; %80, 280, % here things slowly start to converge
        ts5 = '20200522_224032_1_7381'; %50, 330
        ts6='20200523_044512_1_7381';%50, 380
        ts7 = '20200523_044512_2_5536'; %30, 410
        ts8 = '20200523_135237_1_5536'; %50, 460
        ts9 = '20200523_200855_1_3691'; %100, 560
        ts10 = '20200523_200855_2_2214';%30 , 590
        ts11 = '20200524_123605_1_2214';%10, 600
        ts12 = '20200524_140951_1_2214';%20, 620
        ts13 = '20200524_163232_1_2214';%20, 640
        ts14 = '20200524_185508_1_2214';%40, 680
        ts15 = '20200524_233515_1_2214';%100, 780
        ts16 = '20200525_114752_1_1329'; %20, 800
        ts17 = '20200603_044203_1_2214'; %200, changed to 500 walkers. %%cumulativefrom17 200
        ts18 = '20200603_044203_2_1845';%70 %%cumulativefrom17 270
        ts19 = '20200604_073858_1_1845'; %200 %%cumulativefrom17 470
        ts20 = '20200604_073858_2_1476'; % 50  %%cumulativefrom17 520
        ts21 = '20200605_045501_1_1476';% 200  %%cumulativefrom17 720
        ts22 = '20200605_045501_2_1181';% 200 %%cumulativefrom17 920
        ts23 = '20200605_045501_3_1033';% 200 %%cumulativefrom17 1120
        ts24 = '20200605_045501_4_886';% 200 %%cumulativefrom17 1320
        ts25 = '20200605_045501_5_738';% 200 %%cumulativefrom17 1520
        ts26 = '20200605_045501_6_591';% 200  %%cumulativefrom17 1590 + 130
        ts27 = '20200605_045501_7_443';% 200  %%cumulativefrom17 1790 + 130
        ts28 = '20200609_190458_1_1033'; %50  %%cumulativefrom17 1840 + 130
        ts29 = '20200610_000545_1_1033'; %400  %%cumulativefrom17 2240 + 130
        ts30 = '20200610_000545_2_738'; %280  %%cumulativefrom17 2520 + 130
        % v6 starts here. 
        ts31 = '20200612_214710_1_1476'; % 30   %%cumulativefrom17 2550 + 130, windows update after 3 iters
        ts32 = '20200614_005043_1_1476';%100  %%cumulativefrom17 2650 + 130
        ts33 = '20200614_005043_2_2214';%100  %%cumulativefrom17 2750 + 130
        ts34 = '20200614_005043_3_2953';%100  %%cumulativefrom17 2850 + 130
        ts35 = '20200614_005043_4_3691';%100  %%cumulativefrom17 2950 + 130
        ts36 = '20200614_005043_5_2953';%50  %%cumulativefrom17 3000 + 130
        ts37 ='20200615_233435_1_2953';%220  %%cumulativefrom17 3220 + 130
        ts38 = '20200616_183409_1_2953';%400  %%cumulativefrom17 3620 + 130
        ts39 = '20200616_183409_2_2214';%390  %%cumulativefrom17 4010 + 130
        ts40 = '20200619_213706_1_2214';%370  %%cumulativefrom17 4380 + 130
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
        tstamp = ...
            {ts17 ts18 ts19 ts20 ts21 ts22 ts23 ts24 ts25 ts26...
             ts27 ts28 ts29 ts30 ts31 ts32 ts33 ts34 ts35 ts36...
             ts37 ts38 ts39 ts40 ts41 ts42 ts43 ts44 ts45 ts46 ...
             ts47 ts48 ts49 ts50 ts51 ts52 ts53 ts54 ts55 ts56 ...
             ts57 ts58 ts59};
        nIterID = ...
            {1:20 1:7 1:20 1:5 1:20 1:20 1:20 1:20 1:20 1:20 ...
             1:20 1:5 1:40 1:28 1:3 1:10 1:10 1:10 1:10 1:5 ...
             1:22 1:40 1:39 1:40 1:2 1:10 1:36 1:100 1:40 1:40 ...
             1:40 1:29 1:40 1:40 1:40 1:40 1:100 1:21 1:100 1:60 ...
             1:100 1:100 1:3};
         % Just plot v6. if you want to plot everything from ts 17, use the
         % list above. 
        tstamp = ...
            {ts31 ts32 ts33 ts34 ts35 ts36...
             ts37 ts38 ts39 ts40 ts41 ts42 ts43 ts44 ts45 ts46 ...
             ts47 ts48 ts49 ts50 ts51 ts52 ts53 ts54 ts55 ts56 ...
             ts57 ts58 ts59};
        nIterID = ...
            {1:3 1:10 1:10 1:10 1:10 1:5 ...
             1:22 1:40 1:39 1:40 1:2 1:10 1:36 1:100 1:40 1:40 ...
             1:40 1:29 1:40 1:40 1:40 1:40 1:100 1:21 1:100 1:60 ...
             1:100 1:100 1:3};        
        load([projdir '\simdata_' tstamp{end} '\full_variable_set_'...
            tstamp{end} '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info',  'ri');
end
tsToSave = tstamp{end};
mai.masterVector
marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray =  marray_full;%cat(3, marray1, marray_full(:,:,1:end));

clear marray_full
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

mcmc_plot(marray(:, 1:end,1:end), parnames(:),...
    'savematlabfig', false,'plotDistribution', false,...
    'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave);    

%  [activeNames(estParamsIX,1) mat2cell(log(cell2mat(activeNames(estParamsIX,3))), ones(14,1), ones(1, 2))]
%
% ans =
%
%   14×3 cell array
%
%     {'TX_elong_glob'                  }    {[  0]}    {[ 5]}
%     {'TL_elong_glob'                  }    {[  0]}    {[ 6]}
%     {'TXTL_PTET_RNAPbound_Kd'         }    {[  0]}    {[22]}
%     {'TXTL_PLAC_RNAPbound_Kd'         }    {[  5]}    {[17]}
%     {'TXTL_RNAPBOUND_TERMINATION_RATE'}    {[  0]}    {[12]}
%     {'TXTL_RIBOBOUND_TERMINATION_RATE'}    {[  0]}    {[12]}
%     {'RNAP'                           }    {[  1]}    {[15]}
%     {'Ribo'                           }    {[  1]}    {[15]}
%     {'TXTL_PLAS_RNAPbound_Kd'         }    {[ 25]}    {[40]}
%     {'TXTL_PLAS_TFBIND_Kd'            }    {[  0]}    {[10]}
%     {'TXTL_PLAS_TFRNAPbound_Kd'       }    {[  0]}    {[15]}

% if plotflag
% % % %
% % % % close all
% % % % % %    close all
% % % % %%     % Plot trace and corner (posterior distribution) plots
%     mcmc_plot(marray(:, 1:end,1:4:end), parnames(:),...
%         'savematlabfig', figsave, 'savejpeg', jpgsave,...
%         'projdir', projdir, 'tstamp', tsToSave,...
%         'extrafignamestring', 'Without_transient');
% % % % %%
%
%%
mcmc_plot(marray(:, 1:20:end,1:end), parnames(:),...
    'savematlabfig', false, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave);
%%
mcmc_plot(marray(:, 1:end,end-600:100:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'plotChains', false,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'Burned_in_oversampled');
%%
mcmc_plot(marray(:, 1:end,end), parnames(:),...
    'ess', 100, 'scatter', false, 'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'Burned_in_ess100');
% %%
% mcmc_plot(marray(:, 1:20:end,15:end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'fewChains_all_iters_');

%%
mcmc_plot(marray(:, 1:50:end,1:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'fewChains_all_iters');

% % % % %
%%
%
%%


% marray_thin = marray(:, :, end-600:60:end);
% marray_flat = marray_thin(:,:)';
% 
% 
% parslices = [...
%     {'TX_{cat}'}        {[0 6]}
%     {'TL_{cat}'}        {[0 6]}
%     {'\tau_{atp}'}      {[9.6 9.9]}
%     {'pol_{Kd, tet}'}   {[17 23]}
%     {'rep_{Kd, tet}'}   {[-4 1]}
%     {'aTc_{Kd}'}        {[-10 5]}
%     {'pol_{Kd, lac}'}   {[16 23]}
%     {'pol_{term}'}      {[2 10.5]}
%     {'Ribo_{term}'}     {[2 10]}
%     {'pol'}             {[4 8]}
%     {'Ribo'}            {[2 12]}
%     {'3OC12_{Kd}'}      {[5 16]}
%     {'pol_{Kd,las}'}    {[25 40]}
%     {'plas_{tf, Kd}'}   {[4 8]}
%     {'plas-pol_{tf, Kd}'}    {[3 10]}];
% % pol kd is cartesian with pol term, rnaseKd, rnase cat. BUT NOT POL. 
% % cartesian with  rnaseKd, rnase cat and POL
% cutted_marray = mcmc_cut(marray_flat, 1:15, flipud(cell2mat(parslices(:, 2))'));
% size(cutted_marray)

%
% mcmc_plot(cutted_marray(:,:), parnames(:))




%
% parIDs = [4 5 11];
% %     parRanges(parIDs, :) = [-1 0 ;
% %         -4 0;
% %         11 15]
%     parRanges(parIDs, :) = [-1 0 ;
%         -3 -1;
%         11.5 14.5]
%
%     marray_cut = mcmc_cut(marray(:, 1:end,(end - 40):end), parIDs, flipud((parRanges(parIDs, :))'));
%
%     mcut0 = marray_cut(parIDs, :,:);
%     mcmc_3D(mcut0(:,:)', parnames(parIDs), '3 params');
%
%
%     mcmc_plot(marray_cut, parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', false,...
%     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
%     mcut = marray_cut(parIDs, 1:end,ceil(end/4):end);
%     mcmc_3D(mcut(:,:)', parnames(parIDs), '3 params')
%



% We can literally pick any point in this cartesian product.
%     parRanges(parIDs, :) = [-1.16, -0.92 ;
%         -4 -0.5;
%         -15 -8.5]


%     parRanges([1, 2, 3], :) = [15.59 15.73 ;
%         -0.2966 0.1031;
%         8.385 8.622];
%     marray_cut = mcmc_cut(marray, [6,7,9], flipud((parRanges([6 7 9], :))'));
%
%     mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', false,...
%     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
%     mcut = marray_cut([6 7 9], 1:end,ceil(end/4):end);
%     mcmc_3D(mcut(:,:)', parnames([6 7 9]), 'RNA deg covariation')

% Here, we try to plot the trajectories for different sets of parameters.
% I will basically cut the parameters to different subsets, and plot the
% results to check if the lasR model is even working. if it is, then it
% will be clear that it just needs more simulation / fixing of the F rates
% and an estimation of the Kds.
%



% % % %
% % % %     paramIndices = (1:12);
% % % %     parRanges(paramIndices, :) = [...
% % % %         10.5 20.5; %1 pol_{Kd,tet}
% % % %         7 13; %2 pol_{Kd,lac}
% % % %         4 6.1; %3 pol % 5.9 in mcmc_info_ZSIFFL_training_full (from mtet phase 1, 2 etc)
% % % %         5.7 12; %4 Ribo % 5.9 in mcmc_info_ZSIFFL_training_full (from mtet phase 1, 2 etc)
% % % %        -20 20 ; %5 3OC12_{Kd} % correlated with plas_{tf, F}
% % % %         0 2; %6 3OC12_{F} % just gaussian.
% % % %         -20 20; %7 pol_{Kd,las}
% % % %         -2 2.2; %8 pol_{F,las}
% % % %         -20.5 20;%9 plas_{tf, Kd} <
% % % %         -20 20; % 10 -- plas-pol_{tf, Kd} <
% % % %         -2 5;% plas-pol_{tf, F}
% % % %         -2.5 1.5;% plas_{tf, F}
% % % %         ];
% % % %     marray_cut = mcmc_cut(marray, paramIndices, flipud((parRanges(paramIndices, :))'));
% % % %     size(marray_cut)

% actually i think it is better to just set te parameters using gaussian balls,
% and fix as many as possible, and jsut explore the Plas and OC12 Kds.
%
%     % mean and sd parameter.
%     gaussianMeanSD = ...
%         [16.57 0.1; ...%
%         10.5 0.1; ...% pol_{Kd,lac}
%         5.9 0.1; ...% pol
%         5.9 0.1; ...% Ribo
%         12 1; ...% 3OC12_{Kd}
%         0 1 ; ...% 3OC12_{F}
%         35 4; ...% pol_{Kd,las}
%         3 2; ...% pol_{F,las}
%         5 5; ...% plas_{tf, Kd} <
%         8 5; ...% plas-pol_{tf, Kd} <
%         0 2; ...% plas-pol_{tf, F}
%         0 4; ...%plas_{tf, F}
%         ];
%
%    marray_gauss =  cell2mat(arrayfun(@(mu, sig) mu + sig*randn(1, 100),...
%        gaussianMeanSD(:, 1), gaussianMeanSD(:, 2), ...
%        'UniformOutput', false));
%
%


% Plot trajectories.

% rebuild the master vector array, either via mcmc_cut or just using
% all estimated points.
%%
%
%         mvarray = masterVecArray(marray_gauss, mai);
%          mvarray = masterVecArray(marray_cut, mai);
mvarray = masterVecArray(marray, mai);
%     clear marray
%
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
    titls_array
    samplePoints = ceil(size(mvarray, 3) * [1]);
    
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
    
    if finafigmode == false
        mcmc_trajectories(currmi.emo, currdi, currmi, marrayOrd,...
            titls_array, {},...
            'SimMode', 'meanstd', 'separateExpSim', true,...
            'savematlabfig', figsave, 'savejpeg', jpgsave,...
            'projdir', projdir, 'tstamp', tsToSave,...
            'extrafignamestring', ...
            [extrastring num2str(miID) ' ' num2str(msID)]);%,'collateDoses', true,
    end
    
end
%
% end
marrayOrd(:,1:5,end)
[(1:42)' mvarray(:,1:5,end)]
% clear marrayOrd
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]

%% now make the collated trajectory figure and the endpoints figure for the paper.
% also, replot the markov chain trace plots and the posterior distribution
% plots here with proper font sizes.

% resimulate, but with just the matured degfp.
ms = {{'protein deGFP*'}};
for miID = 1:5
    currmi = mi(miID);
    currdi = data_info(currmi.dataToMapTo);
    tv = currdi.timeVector;
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
    dose  = currmi.dosedVals'
    vi = get(currmi.emo, 'ValueInfo')
    for i = 1:length(vi)
        vi(i)
    end
    currmi.dosedNames
    currmi.dosedVals
    
    dose
    
    [da{miID}, idxnotused{miID}] = ...
        simulatecurves(currmi.emo,marrayOrd(:,:)', 50, dose, tv, ms);
end


%% here


ylims = {[7.5, 2, 2, 5, 7.5];
    [7.5, 2, 2, 5, 7.5];
    [15 5 5 10 15]};


lengthToPlotArray = [24, 31, 61]; % 180 min = 24
for outercount = 1:length(lengthToPlotArray)
    lengthToPlot = lengthToPlotArray(outercount)
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
    
    
    
    %     clear marray
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
    

%     print([saveFinalFigs tsToSave 'fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-depsc')
%     print([saveFinalFigs tsToSave 'fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-djpeg')
%     print([saveFinalFigs tsToSave '_lowres_traj_' num2str((lengthToPlot-1)*8) 'min'],'-dpng', '-r200')
%     print([saveFinalFigs tsToSave '_medres_traj_' num2str((lengthToPlot-1)*8) 'min'],'-dpng', '-r300')    
%     print([saveFinalFigs tsToSave '_hires_traj_' num2str((lengthToPlot-1)*8) 'min'],'-dpng', '-r400')
%      
%     print([saveFinalFigs tsToSave '_hires_traj_' num2str((lengthToPlot-1)*8) 'min'],'-dpdf')
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
    %
    
    % Now generate the plots
    % close all
    
    %
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
        ax.XLim = xLims(count, :);
        hold on
        ax = gca
        ax.XTick = fliplr(doseArray(count, 1:2:end));
        ax.XTickLabel = fliplr(xTickLabels(count, :));
        ax.FontSize = 16
        
        %     title(titleArray{count})
        ylabel('deGFP, uM');
        xlabel(xLab{count});
        
    end
%     print([saveFinalFigs tsToSave 'fitting_end_' ...
%         num2str((lengthToPlot-1)*8) 'min'],'-depsc')
%     print([saveFinalFigs tsToSave 'fitting_end_' ...
%         num2str((lengthToPlot-1)*8) 'min'],'-djpeg')
%     print([saveFinalFigs tsToSave 'fitting_end_' ...
%         num2str((lengthToPlot-1)*8) 'min'],'-dpng')
%     print([saveFinalFigs tsToSave '_lowres_end_' ...
%         num2str((lengthToPlot-1)*8) 'min'],'-dpng', '-r200')
%     print([saveFinalFigs tsToSave '_medres_end_' ...
%         num2str((lengthToPlot-1)*8) 'min'],'-dpng', '-r300')
%     print([saveFinalFigs tsToSave '_hires_end_' ...
%         num2str((lengthToPlot-1)*8) 'min'],'-dpng', '-r400')
%     print([saveFinalFigs tsToSave '_hires_end_' ...
%         num2str((lengthToPlot-1)*8) 'min'],'-dpdf')
    
end
%
