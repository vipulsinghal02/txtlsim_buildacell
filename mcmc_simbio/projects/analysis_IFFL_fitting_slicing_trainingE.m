
clear all
close all
saveFinalFigs = '/Users/vipulsinghal/Dropbox/Documents/a_Journal_Papers/Drafts/txtl_bmc_bioinformatics/figs/Jul20_2019/'

finafigmode = true
% Set working directory to the txtlsim toolbox directory.
projdir = '/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/mcmc_simbio/projects/proj_ZSIFFL_trainingE'
addpath(projdir)
sls = regexp(projdir, '/', 'split');
extrastring = sls{end};
jpgsave = false;
figsave = false;

% Load model, mcmc_info, and data_info.
% construct simbiology model object(s)
mtet = model_txtl_ptetdeGFP_pLactetR_aTc;
mlas = model_txtl_pLacLasR_pLasdeGFP;
mlac = model_txtl_pLacdeGFP;
% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_training_fullE(mtet, mlac, mlas);
di = data_ZSIFFL;
mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;
% plot data from existing simulations.
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
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info',  'ri');
end
tsToSave = tstamp{end};
mai.masterVector
marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray = marray_full(:,:,1:end);
clear marray_full
parnames = [...
    {'pol_{Kd, tet}'}
    {'pol_{Kd, lac}'}
    {'pol_{term}'}
    {'Ribo_{term}'}
    {'pol'}
    {'Ribo'}
    {'pol_{Kd,las}'}
    {'plas_{tf, Kd}'}
    {'plas-pol_{tf, Kd}'}    ];

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

%%
%% DIST A
mcmc_plot(marray(:, 1:end,1:10:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'with_transient');

%%

close all
mcmc_plot(marray(:, 1:end,(end - 60):2:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'Burned_in');
% %% DIST B
% parRanges = [17 25
%     17 25
%     2 10.5
%     2 10
%     2.5 10
%     2.5 12
%     25 40
%     0 10
%     0 12]
% parIDs = 1:9;
% marray_cut = mcmc_cut(marray(:, 1:end,(end - 40):end),...
%     parIDs,...
%     flipud((parRanges(parIDs, :))'));
% mcmc_plot(marray_cut, parnames(:));
% %% DIST C
% parRanges = [17 25
%     17 25
%     4.5 8
%     2 10
%     2.5 10
%     2.5 12
%     25 40
%     0 10
%     0 12]
% parIDs = 1:9;
% marray_cut = mcmc_cut(marray(:, 1:end,(end - 40):end),...
%     parIDs,...
%     flipud((parRanges(parIDs, :))'));
% mcmc_plot(marray_cut, parnames(:));
% % % % % %
% %% 
% parRanges = [17 25
%     17 25
%     4.5 8
%     2 10
%     4 8
%     2.5 12
%     25 40
%     4 9
%     3 12]
% parIDs = 1:9;
% marray_cut = mcmc_cut(marray(:, 1:end,(end - 100):end),...
%     parIDs,...
%     flipud((parRanges(parIDs, :))'));
% mcmc_plot(marray_cut, parnames(:));
% %%
% 
% parRanges = [17 25
%     17 25
%     4.5 8
%     2 10
%     4 6.6
%     2.5 12
%     25 40
%     4 9
%     3 12]
% parIDs = 1:9;
% marray_cut = mcmc_cut(marray(:, 1:end,(end - 100):end),...
%     parIDs,...
%     flipud((parRanges(parIDs, :))'));
% % mcmc_plot(marray_cut, parnames(:));
% mcmc_plot(marray_cut, parnames(:), 'ess', 300, 'scatter', false);
% 
% %% and with ribo term truncated at the top: 
% close all
% parRanges = [17 25
%     17 25
%     4.5 8
%     2 6.8
%     4 6.6
%     2.5 12
%     25 40
%     4 9
%     3 12]
% parIDs = 1:9;
% marray_cut = mcmc_cut(marray(:, 1:end,(end - 100):end),...
%     parIDs,...
%     flipud((parRanges(parIDs, :))'));
% 
% mcmc_plot(marray_cut, parnames(:));
% mcmc_plot(marray_cut, parnames(:), 'ess', 300, 'scatter', false);
% 
% %%  DIST D
% %next we tighern all the broad distributions step by step: 
% % first riboterm around its peak
% % ribo 6 to 12
% % pol kd las 29.5 to 40
% close all
% parRanges = [17 25
%     17 25
%     4.5 8
%     5.1 5.7
%     4 6.6
%     6 12
%     29.5 40
%     4 9
%     3 12]
% parIDs = 1:9;
% marray_cut = mcmc_cut(marray(:, 1:end,(end - 100):end),...
%     parIDs,...
%     flipud((parRanges(parIDs, :))'));
% 
% mcmc_plot(marray_cut, parnames(:));
% % mcmc_plot(marray_cut, parnames(:), 'ess', 300, 'scatter', false);
% %% DIST E finally, remove the covariation and look at the resulting spread in the 
% % parameters
% close all
% parRanges = [19.5 20.1
%     17 25
%     4.5 8
%     5.1 5.7
%     4.8 6.4
%     6 12
%     29.5 40
%     4 9
%     3 12]
% parIDs = 1:9;
% marray_cut = mcmc_cut(marray(:, 1:end,(end - 100):end),...
%     parIDs,...
%     flipud((parRanges(parIDs, :))'));
% 
% mcmc_plot(marray_cut, parnames(:));
% % mcmc_plot(marray_cut, parnames(:), 'ess', 300, 'scatter', false);

%%
% DIST F from DIST B, but with covariation removed. 

parRanges = [19.5 20.1
    17 25
    2 10.5
    2 10
    2.5 10
    2.5 12
    25 40
    0 10
    0 12]
parIDs = 1:9;
marray_cut = mcmc_cut(marray(:, 1:end,(end - 100):end),...
    parIDs,...
    flipud((parRanges(parIDs, :))'));
mcmc_plot(marray_cut, parnames(:));


%% DIST G: use the pol restriction from 1.44 to 1.62 used in the 
% vnprl <-- ACTUALLY THIS WIL NOT WORK SINCE tHERE IS NO DENSITY ThERE> 
% SEEMS LIKE A PROBLEM. WHY IS THERE NO DENSITY WHERE THERE SHOULD BE
% ACCORDING TO VNPRL. 

%% 
% Next, I just 
% 1. list out the full parameter values (from the master vector) for all
% the parameters, with the points draws from this distribution. 
% 2. plot the trajectories of the IFFL with these restrictions. 
% 
% FInally, these will be used to generate realistic resutls for the IFFL
% example (tutorial 3), to show how the component config files are
% populated (tutorial 8, not yet)


marray_flat = marray_cut(:,:)';
mvarray = repmat(mai.masterVector, 1, size(marray_flat', 2));
size(marray_flat)
%

estParamsIx = setdiff((1:length(mai.masterVector))', mai.fixedParams);
mvarray(estParamsIx, :) = marray_flat';
size(mvarray)
%%



% Plot trajectories.

% rebuild the master vector array, either via mcmc_cut or just using
% all estimated points.

%
%         mvarray = masterVecArray(marray_gauss, mai);
%          mvarray = masterVecArray(marray_cut, mai);
% mvarray = masterVecArray(marray, mai);
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
    samplePoints = ceil(size(mvarray, 3) * [.95, 1]);
    
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
%%
% end
marrayOrd(:,1:5,end)
[(1:42)' mvarray(:,1:5,end)]
% clear marrayOrd
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]

% now make the collated trajectory figure and the endpoints figure for the paper.
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


% here
close all

ylims = {[7.5, 2, 2, 5, 7.5];
    [7.5, 2, 2, 5, 7.5];
    [15 5 5 10 15]};


lengthToPlotArray = [24, 31, 61]; % 180 min = 24
for outercount = 1:length(lengthToPlotArray)
    lengthToPlot = lengthToPlotArray(outercount)
    
    figure
    ss = get(0, 'screensize');
    set(gcf, 'Position', [ss(3)*(1-1/1.3) ss(4)*(1-1/1.3) ss(3)/3.5 ss(4)/1.1]);
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
    
%     print([saveFinalFigs 'fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-depsc')
%     print([saveFinalFigs 'fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-djpeg')
%     
%     print([saveFinalFigs 'fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-dpng')
    
end


%%

% 
% [activeNames(:,1) mat2cell(exp(mvarray(:,1:8)), ones(42,1), ones(1, 8))]
% 
% ans =
% 
%   42×9 cell array
% 
%     {'TX_elong_glob'    }    {[    22.669]}    {[    22.669]}    {[    22.669]}    {[    22.669]}    {[    22.669]}    {[    22.669]}    {[    22.669]}    {[    22.669]}
%     {'TL_elong_glob'    }    {[    31.062]}    {[    31.062]}    {[    31.062]}    {[    31.062]}    {[    31.062]}    {[    31.062]}    {[    31.062]}    {[    31.062]}
%     {'AGTPdeg_time'     }    {[     23156]}    {[     23156]}    {[     23156]}    {[     23156]}    {[     23156]}    {[     23156]}    {[     23156]}    {[     23156]}
%     {'AGTPreg_ON'       }    {[      0.02]}    {[      0.02]}    {[      0.02]}    {[      0.02]}    {[      0.02]}    {[      0.02]}    {[      0.02]}    {[      0.02]}
%     {'AGTPdeg_rate'     }    {[ 5.616e-05]}    {[ 5.616e-05]}    {[ 5.616e-05]}    {[ 5.616e-05]}    {[ 5.616e-05]}    {[ 5.616e-05]}    {[ 5.616e-05]}    {[ 5.616e-05]}
%     {'TXTL_UTR_UTR1_Kd' }    {[    1.0557]}    {[    1.0557]}    {[    1.0557]}    {[    1.0557]}    {[    1.0557]}    {[    1.0557]}    {[    1.0557]}    {[    1.0557]}
%     {'TXTL_PTET_RNAPb?'}    {[3.2394e+08]}    {[3.4631e+08]}    {[5.2969e+08]}    {[5.3155e+08]}    {[3.4911e+08]}    {[3.7055e+08]}    {[3.0805e+08]}    {[3.9848e+08]}
%     {'TXTL_PTET_RNAPb?'}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}
%     {'TXTL_NTP_RNAP_1?'}    {[    19.028]}    {[    19.028]}    {[    19.028]}    {[    19.028]}    {[    19.028]}    {[    19.028]}    {[    19.028]}    {[    19.028]}
%     {'TXTL_NTP_RNAP_2?'}    {[ 1.199e+06]}    {[ 1.199e+06]}    {[ 1.199e+06]}    {[ 1.199e+06]}    {[ 1.199e+06]}    {[ 1.199e+06]}    {[ 1.199e+06]}    {[ 1.199e+06]}
%     {'TXTL_PTET_seque?'}    {[   0.60653]}    {[   0.60653]}    {[   0.60653]}    {[   0.60653]}    {[   0.60653]}    {[   0.60653]}    {[   0.60653]}    {[   0.60653]}
%     {'TXTL_PTET_seque?'}    {[     3.721]}    {[     3.721]}    {[     3.721]}    {[     3.721]}    {[     3.721]}    {[     3.721]}    {[     3.721]}    {[     3.721]}
%     {'TL_AA_Kd'         }    {[    703.87]}    {[    703.87]}    {[    703.87]}    {[    703.87]}    {[    703.87]}    {[    703.87]}    {[    703.87]}    {[    703.87]}
%     {'TL_AGTP_Kd'       }    {[2.0007e+06]}    {[2.0007e+06]}    {[2.0007e+06]}    {[2.0007e+06]}    {[2.0007e+06]}    {[2.0007e+06]}    {[2.0007e+06]}    {[2.0007e+06]}
%     {'TXTL_RNAdeg_Kd'   }    {[6.1681e+06]}    {[6.1681e+06]}    {[6.1681e+06]}    {[6.1681e+06]}    {[6.1681e+06]}    {[6.1681e+06]}    {[6.1681e+06]}    {[6.1681e+06]}
%     {'TXTL_INDUCER_TE?'}    {[   0.13534]}    {[   0.13534]}    {[   0.13534]}    {[   0.13534]}    {[   0.13534]}    {[   0.13534]}    {[   0.13534]}    {[   0.13534]}
%     {'TXTL_INDUCER_TE?'}    {[    4.8404]}    {[    4.8404]}    {[    4.8404]}    {[    4.8404]}    {[    4.8404]}    {[    4.8404]}    {[    4.8404]}    {[    4.8404]}
%     {'TXTL_DIMER_tetR?'}    {[  4.54e-05]}    {[  4.54e-05]}    {[  4.54e-05]}    {[  4.54e-05]}    {[  4.54e-05]}    {[  4.54e-05]}    {[  4.54e-05]}    {[  4.54e-05]}
%     {'TXTL_DIMER_tetR_F'}    {[    4.2503]}    {[    4.2503]}    {[    4.2503]}    {[    4.2503]}    {[    4.2503]}    {[    4.2503]}    {[    4.2503]}    {[    4.2503]}
%     {'TXTL_UTR_UTR1_F'  }    {[   0.81873]}    {[   0.81873]}    {[   0.81873]}    {[   0.81873]}    {[   0.81873]}    {[   0.81873]}    {[   0.81873]}    {[   0.81873]}
%     {'TXTL_PLAC_RNAPb?'}    {[2.5133e+08]}    {[1.9969e+08]}    {[5.4477e+08]}    {[4.3346e+08]}    {[2.9277e+08]}    {[3.8904e+08]}    {[2.6869e+08]}    {[2.6059e+08]}
%     {'TXTL_PLAC_RNAPb?'}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}    {[    4.4817]}
%     {'TXTL_RNAPBOUND_?'}    {[    238.48]}    {[    4991.1]}    {[    1702.3]}    {[    284.84]}    {[    8476.3]}    {[     981.5]}    {[    12.086]}    {[    217.52]}
%     {'TXTL_NTP_RNAP_1_F'}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}
%     {'TXTL_NTP_RNAP_2_F'}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}
%     {'TL_AA_F'          }    {[   0.74082]}    {[   0.74082]}    {[   0.74082]}    {[   0.74082]}    {[   0.74082]}    {[   0.74082]}    {[   0.74082]}    {[   0.74082]}
%     {'TL_AGTP_F'        }    {[   0.30119]}    {[   0.30119]}    {[   0.30119]}    {[   0.30119]}    {[   0.30119]}    {[   0.30119]}    {[   0.30119]}    {[   0.30119]}
%     {'TXTL_RIBOBOUND_?'}    {[    501.46]}    {[      1528]}    {[    192.79]}    {[    255.95]}    {[    162.69]}    {[    81.165]}    {[    1554.2]}    {[    499.89]}
%     {'TXTL_RNAdeg_F'    }    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}
%     {'TXTL_RNAdeg_kc'   }    {[   0.79844]}    {[   0.79844]}    {[   0.79844]}    {[   0.79844]}    {[   0.79844]}    {[   0.79844]}    {[   0.79844]}    {[   0.79844]}
%     {'RNAP'             }    {[    157.02]}    {[    160.85]}    {[    423.81]}    {[    383.77]}    {[     273.5]}    {[    330.85]}    {[     210.8]}    {[     244.7]}
%     {'RNase'            }    {[    9897.1]}    {[    9897.1]}    {[    9897.1]}    {[    9897.1]}    {[    9897.1]}    {[    9897.1]}    {[    9897.1]}    {[    9897.1]}
%     {'Ribo'             }    {[     11660]}    {[    8170.5]}    {[    1340.6]}    {[    5572.7]}    {[    9732.7]}    {[     21461]}    {[    6322.9]}    {[     41499]}
%     {'TXTL_PROT_deGFP?'}    {[ 0.0023001]}    {[ 0.0023001]}    {[ 0.0023001]}    {[ 0.0023001]}    {[ 0.0023001]}    {[ 0.0023001]}    {[ 0.0023001]}    {[ 0.0023001]}
%     {'TXTL_INDUCER_LA?'}    {[4.4241e+05]}    {[4.4241e+05]}    {[4.4241e+05]}    {[4.4241e+05]}    {[4.4241e+05]}    {[4.4241e+05]}    {[4.4241e+05]}    {[4.4241e+05]}
%     {'TXTL_INDUCER_LA?'}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}
%     {'TXTL_PLAS_RNAPb?'}    {[6.8311e+14]}    {[9.9043e+13]}    {[2.8567e+13]}    {[4.5872e+14]}    {[1.0489e+12]}    {[3.7514e+13]}    {[8.1248e+15]}    {[7.5184e+14]}
%     {'TXTL_PLAS_RNAPb?'}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}
%     {'TXTL_PLAS_TFBIN?'}    {[    142.43]}    {[    73.651]}    {[    345.12]}    {[    1234.3]}    {[    669.89]}    {[    28.171]}    {[    609.41]}    {[    760.69]}
%     {'TXTL_PLAS_TFRNA?'}    {[    2808.8]}    {[     37613]}    {[     14653]}    {[      1060]}    {[    4084.6]}    {[    133.54]}    {[    3238.9]}    {[     172.1]}
%     {'TXTL_PLAS_TFRNA?'}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}
%     {'TXTL_PLAS_TFBIN?'}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}    {[         1]}
% 
