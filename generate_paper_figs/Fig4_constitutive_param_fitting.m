% File: Fig4_constitutive_param_fitting.m
%
% this file is a clone of analysis_vnprl_F2_v2.m. It generates Fig 4, and
% Supp fig 4. 
%
% Purpose: Plot the VNPRL 2011 + ACS 2014 papers data.
%
% Plot the fits to that data that were performed using mcmc_simbio.
%
% Author:
% Vipul Singhal
%
% The relevant literature for this data is:
%
% 1. Gene Circuit Performance Characterization and Resource Usage in a 
% Cell-Free Breadboard, Siegal-Gaskins, et. al
%
% 2. Coarse-Grained Dynamics of Protein Synthesis in a Cell-Free System
% Karzbrun et. al


% ----------------------------------------------------------------------- %
% Set data and input directories
% ----------------------------------------------------------------------- %

% set flags for whether to save the plots generated (default is false --
% images are large, make sure there is enough RAM / Hard disk space)
jpgsave = false;

% Set the working directory to be the
txtldir = txtl_init; % txtl_init is in the trunk, so this sets everything
mcmc_init;

projdir = [pwd '/mcmc_simbio/projects/proj_vnprl'];

addpath(projdir)

% this dataset is ~5GB. So worth storing it on an external drive. Change
% the directory below accordingly. 
external_drive_where_data_is_stored = '/Volumes/vs_mark1/CALTECH/2020_FEB/proj_vnprl';
if ~isempty(external_drive_where_data_is_stored)
    data_dir = external_drive_where_data_is_stored;
    addpath(data_dir)
else
    data_dir = projdir;
end


% ----------------------------------------------------------------------- %
% Load the models and the data
% ----------------------------------------------------------------------- %

% Load model, mcmc_info, and data_info.
mobj = model_dsg2014_regen;
mcmc_info = mcmc_info_vnprl_F2(mobj);
di = data_VNPRL2011; % change this.

mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

% plot data from existing simulations.

ts1 = '20190218_135635_1_24';
ts2 = '20190218_135635_2_12';
ts3 = '20190219_093714_1_7';
ts4 = '20190219_093714_2_2';
ts5 = '20190220_022918_1_2';
ts6 = '20190220_102006_1_2';
ts7 = '20190221_074748_1_2';
ts8 = '20190223_024333_1_2';

tstamp = {ts1 ts2 ts3 ts4 ts5 ts6 ts7 ts8}; %
nIterID = {1:25 1:25 1:20 1:20 1:16 1:31 1:80 1:20}; %
load([data_dir '/simdata_' ts8 '/full_variable_set_' ts8 '.mat'], ...
    'mi',...
    'mcmc_info', 'data_info', 'mai', 'ri');

tsToSave = ts8;
mai.masterVector

marray_full = mcmc_get_walkers(tstamp,nIterID, data_dir);
marray = marray_full(:,:,1:end);
clear marray_full
parnames = ...
    [{'TX_{cat}'    }
    {'\tau_{atp}'   }
    {'\delta_{atp}' }
    {'pol_{Kd}'     }
    {'pol_{term}'   }
    {'RNAse_{Kd}'   }
    {'RNAse_{cat}'  }
    {'pol'          }
    {'RNase'        }
    {'TL_{cat}'     }
    {'Ribo_{Kd}'    }
    {'Ribo_{term}'  }
    {'Ribo'         }];



% ----------------------------------------------------------------------- %
% Plot the parameter distributions, and Markov chains.
% ----------------------------------------------------------------------- %

% Uncomment thiese to plot other variations of the corner plot and the
% trace plots:

% mcmc_plot(marray(:, 1:10:end,1:1000:end), parnames(:),...
%     'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', tsToSave, ...
%     'extrafignamestring','WithTransient_jul31');

% mcmc_plot(marray(:, 1:end,(end-10000):4800:end), parnames(:),...
%     'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', tsToSave, ...
%     'extrafignamestring', 'BurnedIn_jul31');
%
mcmc_plot(marray(:, 1:end,(end-10000):3000:end), parnames(:),...
    'scatter', false,'ess', 100, 'plotChains', false,...
    'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, ...
    'extrafignamestring', 'BurnedIn_kde');

%
% ----------------------------------------------------------------------- %
% Plot the trajectories.
% ----------------------------------------------------------------------- %

workingDir = [pwd '/mcmc_simbio/exp_data/public_data/'];

% run merger
%all_data_merger;
load([workingDir 'mergedExperimentFiles.mat'])

colorCodes = {'r','b','g','c','m','k',[1 1 .5],[.7 .5 .2],[0 1 .2],[.35 .8 .8],[.9 0 .4],[1 .2 .2]};

RFUumConvert = 1.723;   % (1.723 a.u. = 1 nM)
RFUumConvertMG = 7.75;  % (7.75 a.u. = 1 nM)


ylims = [1000 50 12000];
lengthToPlotArray = [21, 76, 76];

figure
ss = get(0, 'screensize');
set(gcf, 'Position', [ss(3)*(1-1/1.3) ss(4)*(1-1/1.3) ss(3)/3 ss(4)/1.4]);
miToUse = [1 2 2]; % mi 1 has one plot (RNA) and mi 2 has 2 plots (RNA and protein)
msToPlot = [1 1 2];
timeinterval = [6 8 8];
%
titleArray = {'RNA deg';
    'mRNA expression';
    'deGFP expression'};
legLoc = {'NorthEast','NorthEast','NorthWest'};

legends = [{fliplr({'1000nM', '800nM', '600nM', '200nM', '75nM', '37.5nM'})};
    {{'0.5nM', '2nM', '5nM', '20nM'}};
    {{'0.5nM', '2nM', '5nM', '20nM'}};];
dosesToPlot = {{[1 2 4 5 7 9]}
    {[1 2 3 4]}
    {[1 2 3 4]}};

%     clear marray
yLab = {'RNA, nM', 'RNA, nM', 'deGFP, nM'}

mvarray = masterVecArray(marray(:,:,1:50:end), mai);

samplePoints = ceil(size(mvarray, 3) * [.9, 1]);
    
% resimulate, 
% ms = {{'protein deGFP*'}};
for miID = 1:length(mi)
    
    currmi = mi(miID);
    ms = currmi.measuredSpecies;
    currdi = data_info(currmi.dataToMapTo);
    tv = currdi.timeVector;
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
    dose  = currmi.dosedVals';
    vi = get(currmi.emo, 'ValueInfo');
    for i = 1:length(vi)
        vi(i)
    end
    currmi.dosedNames
    currmi.dosedVals
    
    dose
    
    [da{miID}, idxnotused{miID}] = ...
        simulatecurves(currmi.emo,marrayOrd(:,:)', 50, dose, tv, currmi.measuredSpecies);
end

t_vec_mins=mergedExpFile(1).t_vec/60;

for count = 1:length(miToUse)  %1:length(mi)%1:
    colorz = flipud(parula(max(dosesToPlot{count}{:})+2));
    currmi = mi(miToUse(count));
    currdi = data_info(currmi.dataToMapTo);
    tv = currdi.timeVector;
    lengthToPlot = lengthToPlotArray(count);
    
    % experimental data
    subplot(3, 2, (count-1)*2+1)
    if count == 1 % for the RNA use the data info.
        %
        for dID = dosesToPlot{count}{:}
            plot(tv(1:lengthToPlot)/60,...
                mean(currdi.dataArray(1:lengthToPlot, msToPlot(count), 1, dID),3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            
            hold on
        end
        grid on
        title([titleArray{count} ' (Exp)'])
        ylabel(yLab{count}, 'FontSize', 14)
        axis([0 (lengthToPlot-1)*timeinterval(count), 0 ylims(count)])
        ax = gca;
        ax.FontSize = 14;
        
    elseif count == 2 % for the MG aptamer and the deGFP, use the scripts provided by Zoltan.
        %%%%%  MG KINETICS WITH CORRECT BACKGROUNDS SUBTRACTED
        plotSelect=[ 5, 7, 8, 10];
        mgBackground = mergedExpFile(1).Data_mean(:,1:11,1);
        mgNoBgMean=mergedExpFile(1).Data_mean(:,12:22,1)-mgBackground;
        for i=1:11
            mgNoBgMean(:,i)=mean([ ...
                expFile(1).Data(:,12+i,1)-mgBackground(:,i) ...
                expFile(2).Data(:,12+i,1)-mgBackground(:,i) ...
                expFile(2).Data(:,35+i,1)-mgBackground(:,i) ...
                ],2)/RFUumConvertMG/10;
            
            mgNoBgErr(:,i)=std([ ...
                expFile(1).Data(:,12+i,1)-mgBackground(:,i) ...
                expFile(2).Data(:,12+i,1)-mgBackground(:,i) ...
                expFile(2).Data(:,35+i,1)-mgBackground(:,i) ...
                ],0,2)/sqrt(3)/RFUumConvertMG/10;
        end
%         figure('Name','MG kinetics');
        hold on
        % errorbar(repmat(t_vec_mins,1,length(plotSelect)),mgNoBgMean(:,plotSelect),mgNoBgErr(:,plotSelect));
        for dID = dosesToPlot{count}{:}
            colorzCell{dID} = colorz(dID+2, :);
        end

    
        h = stdshade(t_vec_mins, ...
            mgNoBgMean(:,plotSelect), ...
            mgNoBgErr(:,plotSelect), ...
            0.3, ...
            colorzCell);
        title([titleArray{count} ' (Exp)'])
        %h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,plotSelect+11,1)/RFUumConvertMG,mergedExpFile(1).noBg_std(:,plotSelect+11)/RFUumConvertMG,0.3,colorCodes);
        axis([0 (lengthToPlot-1)*timeinterval(count), ...
            0, ...
            ylims(count)])
        ylabel('RNA, nM')
        grid on
        legendStr = cellstr(num2str(mergedExpFile(1).concentrations(plotSelect)', ...
            'pr-gfp-mg15 %0.2g nM'));
        %legend(h,legendStr,'Location','NorthEast')
        ax = gca;
        ax.FontSize = 14;
    elseif count == 3
        plotSelect=[5, 7, 8, 10];
%         figure('Name','GFP kinetics');
        hold on
        for dID = dosesToPlot{count}{:}
            colorzCell{dID} = colorz(dID+2, :);
        end
        h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,11+plotSelect,2)/RFUumConvert/1.8,mergedExpFile(1).noBg_std(:,11+plotSelect,2)/RFUumConvert/sqrt(3)/1.8,0.3,colorzCell);
        axis([0 (lengthToPlot-1)*timeinterval(count),0,ylims(count)])
        title([titleArray{count} ' (Exp)'])
        xlabel('Time [mins]')
        ylabel('deGFP (uM)')
        grid on
        legendStr = cellstr(num2str(mergedExpFile(1).concentrations(plotSelect)', 'pr-gfp-mg15 %0.2g nM'));
        ax = gca;
        ax.FontSize = 14;
        %legend(h,legendStr,'Location','NorthWest')
    end
    
    
    % simulation results
    subplot(3, 2, (count-1)*2+2)
    for dID = dosesToPlot{count}{:}
        plot(tv(1:lengthToPlot)/60,...
            mean(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID), 3),...
            'LineWidth', 1.5,...
            'Color', colorz(dID+2, :))
        hold on
    end
    grid on
    title([titleArray{count} ' (Fit)'])
    
    legend(legends{count, :}, 'Location', legLoc{count}, 'FontSize', 14)
    legend('boxoff')
    axis([0 (lengthToPlot-1)*timeinterval(count) 0 ylims(count)])
    if count ==length(miToUse)
        xlabel('time, minutes', 'FontSize', 14)
    end
    ax = gca;
    ax.FontSize = 14;
end








% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %