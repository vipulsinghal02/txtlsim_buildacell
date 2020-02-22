


% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_vnprl'];
saveFinalFigs = '/Users/vipulsinghal/Dropbox/Documents/a_Journal_Papers/Drafts/txtl_bmc_bioinformatics/figs/Jul20_2019/'

addpath(projdir)

jpgsave = false;
figsave = false;

% Load model, mcmc_info, and data_info.
mobj = model_dsg2014_regen;
mcmc_info = mcmc_info_vnprl_F2(mobj);
di = data_VNPRL2011;

mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

% plot data from existing simulations.
tsIDtouse = 3;
plotflag = true;
switch tsIDtouse
    case 1
        ts1 = '20190217_154138_1_71';
        ts2 = '20190217_154138_2_48';
        ts3 = '20190217_154138_3_24';
        ts4 = '20190217_154138_4_12';
        tstamp = {ts1 ts2 ts3 ts4};
        nIterID = {1:2 1:2 1:2 1:2};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
    case 2
        ts1 = '20190217_160345_1_48';
        tstamp = {ts1};
        nIterID = {1:9};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');   
    case 3
        ts1 = '20190218_135635_1_24';
        ts2 = '20190218_135635_2_12';
        ts3 = '20190219_093714_1_7';
        ts4 = '20190219_093714_2_2';
        ts5 = '20190220_022918_1_2';
        ts6 = '20190220_102006_1_2';
        ts7 = '20190221_074748_1_2';
        ts8 = '20190223_024333_1_2';
        tstamp = {ts1 ts2 ts3 ts4 ts5 ts6 ts7 ts8};
        nIterID = {1:25 1:25 1:20 1:20 1:16 1:31 1:80 1:20};
        load([projdir '/simdata_' ts8 '/full_variable_set_' ts8 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');  
end
tsToSave = ts8;
mai.masterVector


marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
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
% 
%     {'TX_elong_glob'                  }
%     {'AGTPdeg_time'                   }
%     {'AGTPdeg_rate'                   }
%     {'TXTL_P70_RNAPbound_Kd'          }
%     {'TXTL_RNAPBOUND_TERMINATION_RATE'}
%     {'TXTL_RNAdeg_Kd'                 }
%     {'TXTL_RNAdeg_kc'                 }
%     {'RNAP'                           }
%     {'RNase'                          }
%     {'TL_elong_glob'                  }
%     {'TXTL_UTR_UTR1_Kd'               }
%     {'TXTL_RIBOBOUND_TERMINATION_RATE'}
%     {'Ribo'                           }
%


%% 
%%
if plotflag
    close all
    % Plot trace and corner (posterior distribution) plots
%     mcmc_plot(marray(:, 1:3:end,1:50:end), parnames(:));
%     mcmc_plot(marray(:, 1:5:end,10000:500:end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedInAllWalkers');
% 
%%
    mcmc_plot(marray(:, 1:10:end,1:100:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'WithTransient_jul31');
%%
    mcmc_plot(marray(:, 1:end,(end-12000):600:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn_jul31');
    
%% Stepwise restriction of parameter ranges to get a point that exists in the joint parameter space. 
    % cleaning up the parameters by removing any parameters that are 
    % -2.5, 1.9 ---> % RNAse_{cat} -- 7
    % 3.5, 10.75 ---> % RNase -- 9 
    % because we do not want to be hitting the boundaries of the parameter
    % ranges, there is a biased concentration of points at these
    % boundaries. 
    parRanges = mai.paramRanges;
    parRanges([7,9], :) = [-1.863, 1.9;3.5, 10.75]
    marray_cut = mcmc_cut(marray, [7, 9], flipud((parRanges([7 9], :))'));
    
    mcmc_plot(marray_cut(:, 1:end,10000:500:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');

    % next we restrict the covarying points to small cartesian squares. 
    % RNAse_cat (index 7) ---> -0.2966 0.1031
    % RNase (index 9) ---> 8.298 8.738
    parRanges([7,9], :) = [-0.2966 0.1031; 8.298 8.738];
    marray_cut = mcmc_cut(marray, [7,9], flipud((parRanges([7 9], :))'));
    
    mcmc_plot(marray_cut(:, 1:end,ceil(end/3):10:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
    
    % next we restrict the covarying points to small cartesian squares. 
    % RNAse_Kd (index 6) ---> 15.59 15.73 
    % RNAse_cat (index 7) ---> -0.2966 0.1031
    % RNase (index 9) ---> 8.385 8.622
    parRanges([6,7, 9], :) = [15.59 15.73 ; 
        -0.2966 0.1031;
        8.385 8.622];
    marray_cut = mcmc_cut(marray, [6,7,9], flipud((parRanges([6 7 9], :))'));
    
    mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
    mcut = marray_cut([6 7 9], 1:end,ceil(end/4):end);
    mcmc_3D(mcut(:,:)', parnames([6 7 9]), 'RNA deg covariation')
    
% an attempt at sequential elimination of covariation. The three parameters
% RNase KD, RNase Cat, and RNase covary mutually. We want to find *a*
% volumme hyper cube within the covarying plane. Ideally of maximized mass,
% but otherwise any volume will do. 

% Looks like RNase Kd (15.59 15.73) and RNAse (8.385 8.622) gives a square,
% but the other two covary: the Kd with the cat, and the cat with the
% rnase. 

% as long as I pick values within the first square, it should stay a
% square. So lets try to squarify the other two. 

% KD vs cat
% Kd: 15.65 15.73
% cat: -0.1224 0.02126 

% rnase vs cat
% rnase: 8.386 8.461, 
% cat: -0.039 0.022 --> Kd further restricted to ... nope, the KD vs cat 
% range is cartesian, so no further restriction needed due to the cat being 
% further restricted.
% 

% The intersection of these sets should give us the cube:
% Kd: 15.65 15.73
% cat: -0.039 0.02126 
% rnase: 8.386 8.461 

% lets try this. 
    parRanges([6,7, 9], :) = [15.65 15.73 ; 
        -0.039 0.02126;
        8.386 8.461 ];
    marray_cut = mcmc_cut(marray, [6,7,9], flipud((parRanges([6 7 9], :))'));
    
    mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
    
    % nope, I do not get a cube. Need to think about this. I think it is a
    % tilted hyperplane in space, and just ntersecting a cube with it is
    % not going to give a filled up cube. but even then, its projection
    % should be a bunch of squares, no? need to think about it. this is
    % very interesting. Think about what happens to a tilted plane when you
    % try to intersect a hypercube with it. Also, it is just 3d, can try to
    % just plot it in 3d. oh yeah, lol. 

%     Actually yeah, see the mcmc_3D plot above, it explains this. 


%% FINALLY USED VALUES SECTION %%%%%
% Ok, next up, we are going to restrict most parameters to some ranges, and
% then pick parameter values at random from the set of remaining points. 
% TX cat: 2.3 2.95
    % RNAse_Kd (index 6) ---> 15.59 15.73 
    % RNAse_cat (index 7) ---> -0.2966 0.1031
    % RNase (index 9) ---> 8.385 8.622
    close all
    paramIndices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
    parRanges(paramIndices, :) = [...
        2.3 2.95; %1
        8.6 8.95; %2
        -10.1 -9.7; %3
        13 15; %4
        -1 5 ; %5
        15.59 15.73; %6
        -0.2966 0.1031; %7
        0.7 2.2; %8 pol
        8.385 8.622;%9 RNase
        3.18 3.69; % 10 -- TLcat
        -3 13.5;%RiboKd
        2 3;% Ribo term
        3 4.6 ];% Ribo
    marray_cut = mcmc_cut(marray, paramIndices, flipud((parRanges(paramIndices, :))'));
    
    mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');


 CandidateParams = marray_cut(:,1:100:end,end)
%     figure
%     [C,lags,ESS]=eacorr(marray(:, :,1:end));%10000:end
%     plot(lags,C,'.-',lags([1 end]),[0 0],'k');
%     grid on
%     xlabel('lags')
%     ylabel('autocorrelation');
%     text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',...
%         ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
%     title('Markov Chain Auto Correlation')
%     
%%
%     mvarray = masterVecArray(marray_cut, mai); % instead of marray_cut,
%     lets see what trajectories we get with just marray. 
% this is because marray_cut is, in some sense, biased because we had
% human-in-the-loop eyeballing. 

% CRITICAL: the arrays can be quite large, and eat up all the ram in the
% computer. For this reason, delete arrays as you no longer need them. 
    mvarray = masterVecArray(marray(:,:,1:50:end), mai);
    clear marray
    samplePoints = ceil(size(mvarray, 3) * [.9, 1]);
    %
    marrayOrd = mvarray(mi(1).paramMaps(mi(1).orderingIx),:,samplePoints);
    titls = arrayfun(@num2str, mi(1).dosedVals, 'UniformOutput', false);
    titls_array = cell(length(titls), 2, length(mi(1).measuredSpeciesIndex));
    ms = {'RNA'};
    for i = 1:length(mi(1).measuredSpeciesIndex)
        for j = 1:length(titls)
            titls_array(j, 1, i) = {[ms{i} ', ' titls{j} 'nM initial RNA, Exp data']};
            titls_array(j, 2, i) = {[ms{i} ', ' titls{j} 'nM initial RNA, MCMC samples']};
        end
    end
%     mcmc_trajectories(mi(1).emo, data_info(mi(1).dataToMapTo), mi(1), marrayOrd,...
%         titls_array, {},...
%         'SimMode', 'meanstd', 'separateExpSim', true,...
%         'savematlabfig', figsave, 'savejpeg', jpgsave,...
%         'projdir', projdir, 'tstamp', tsToSave,...
%         'extrafignamestring', 'RNAspike');
    
    marrayOrd = mvarray(mi(2).paramMaps(mi(2).orderingIx),:,samplePoints);
%     clear mvarray
    titls = arrayfun(@num2str, mi(2).dosedVals, 'UniformOutput', false);
    titls_array = cell(length(titls), 2, length(mi(2).measuredSpeciesIndex));
    ms = {'MG aptamer', 'deGFP'};
    for i = 1:length(mi(2).measuredSpeciesIndex)
        for j = 1:length(titls)
            titls_array(j, 1, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, Exp data']};
            titls_array(j, 2, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, MCMC samples']};
        end
    end
    
%     mcmc_trajectories(mi(2).emo, data_info(mi(2).dataToMapTo), mi(2), marrayOrd,...
%         titls_array, {},...
%         'SimMode', 'meanstd', 'separateExpSim', true,...
%         'savematlabfig', figsave, 'savejpeg', jpgsave,...
%         'projdir', projdir, 'tstamp', tsToSave,...
%         'extrafignamestring', 'MGa_deGFP');
end
% marrayOrd(:,1:5,end)
clear marrayOrd
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]

%% remake the trajectory plots so that they are 1. collated and 2. less broad. 
% same type of code as analysis_ZSIFFL_training_fullE.m


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
%%


ylims = [1000 50 12000];

% here, each of the two mi have different length data. 
% the first mi is the spike in rna 
% the second mi is the mrna and gfp expression 
lengthToPlotArray = [21, 76, 76]; 
    
    figure
    ss = get(0, 'screensize');
    set(gcf, 'Position', [ss(3)*(1-1/1.3) ss(4)*(1-1/1.3) ss(3)/3.5 ss(4)/1.6]);
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
    
    
    for count = 1:length(miToUse)  %1:length(mi)%1:
        colorz = flipud(parula(max(dosesToPlot{count}{:})+2));
        currmi = mi(miToUse(count)); 
        currdi = data_info(currmi.dataToMapTo);
        tv = currdi.timeVector;
        lengthToPlot = lengthToPlotArray(count);
        subplot(3, 2, (count-1)*2+1)
        for dID = dosesToPlot{count}{:}
            plot(tv(1:lengthToPlot)/60,...
                mean(currdi.dataArray(1:lengthToPlot, msToPlot(count), 1, dID),3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            
            hold on
        end
        title([titleArray{count} ' (Exp)'])
        ylabel(yLab{count}, 'FontSize', 14)
        axis([0 (lengthToPlot-1)*timeinterval(count), 0 ylims(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
        
        subplot(3, 2, (count-1)*2+2)
        for dID = dosesToPlot{count}{:}
            plot(tv(1:lengthToPlot)/60,...
                mean(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID), 3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            hold on
        end
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
    
%     print([saveFinalFigs 'vnprl_F2_traj_' num2str((lengthToPlot-1)*8) 'min'],'-depsc')
%     print([saveFinalFigs 'vnprl_F2_traj_' num2str((lengthToPlot-1)*8) 'min'],'-djpeg')
%     
%     print([saveFinalFigs 'vnprl_F2_traj_' num2str((lengthToPlot-1)*8) 'min'],'-dpng')
    
% end

%% redo the plot in the previous cell, but with shading of the standard deviation
% the experimental data is actually from Dan Gaskins and Zoltan Tuza's
% paper. 

% The data is stored in the directory
% mcmc_simbio/exp_data/public_data/
% and the code is in the file all_data_compare.m 
% We adapt the code 


    % code we will use is: 
    [linehandle, ptchhandle] =...
        boundedline(tv, summst(:,msi,1, dosei),...
        spreadst(:,msi,1, dosei));
    set(ptchhandle, ...
        'FaceColor', p.SpreadColor,...
        'FaceAlpha', p.FaceAlpha);
    set(linehandle, ...
        'Color', p.LineColor,...
        'LineStyle', p.LineStyle,...
        'LineWidth', p.LineWidth);
    hold on
    
    

ylims = [1000 50 12000];

% here, each of the two mi have different length data. 
% the first mi is the spike in rna 
% the second mi is the mrna and gfp expression 
lengthToPlotArray = [21, 76, 76]; 
    
    figure
    ss = get(0, 'screensize');
    set(gcf, 'Position', [ss(3)*(1-1/1.3) ss(4)*(1-1/1.3) ss(3)/3.5 ss(4)/1.6]);
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
    
    
    for count = 1:length(miToUse)  %1:length(mi)%1:
        colorz = flipud(parula(max(dosesToPlot{count}{:})+2));
        currmi = mi(miToUse(count)); 
        currdi = data_info(currmi.dataToMapTo);
        tv = currdi.timeVector;
        lengthToPlot = lengthToPlotArray(count);
        subplot(3, 2, (count-1)*2+1)
        for dID = dosesToPlot{count}{:}
            
            [linehandle, ptchhandle] = ...
                boundedline(tv(1:lengthToPlot)/60, summst(:,msi,1, dosei),...
                spreadst(:,msi,1, dosei));
            
            plot(tv(1:lengthToPlot)/60, ...
                mean(currdi.dataArray(1:lengthToPlot, msToPlot(count), 1, dID),3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            
            hold on
        end
            [linehandle, ptchhandle] =...
        boundedline(tv, summst(:,msi,1, dosei),...
        spreadst(:,msi,1, dosei));
    
    
        title([titleArray{count} ' (Exp)'])
        ylabel(yLab{count}, 'FontSize', 14)
        axis([0 (lengthToPlot-1)*timeinterval(count), 0 ylims(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
        
        subplot(3, 2, (count-1)*2+2)
        for dID = dosesToPlot{count}{:}
            plot(tv(1:lengthToPlot)/60,...
                mean(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID), 3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            hold on
        end
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



%% find a parameter point that gives good looking fits

% Reaction Strings

%     [{[    6.0913]}    {'TX_elong_glob'                  }  {'[CUTP:AGTP:RNAP:DNA p70--utr1--deGFP] -> [term_RNAP:DNA p70--utr1--deGFP] + [RNA utr1--deGFP]'}
%     {[   24.4062]}    {'TL_elong_glob'                  }   {'[AA:AGTP:Ribo:RNA utr1--deGFP] -> [term_Ribo:RNA utr1--deGFP] + [protein deGFP] + AGMP'}
%     {[5.7776e+03]}    {'AGTPdeg_time'                   }   {'global'}
%     {[    0.0200]}    {'AGTPreg_ON'                     }   {'global'}
%     {[4.8661e-05]}    {'AGTPdeg_rate'                   }   {'AGTP -> AGMP'}
%     {[    0.0023]}    {'TXTL_PROT_deGFP_MATURATION'     }   {'[protein deGFP] -> [protein deGFP*]'}
%     {[  421.4700]}    {'TXTL_UTR_UTR1_Kd'               }   {''}
%     {[    0.8187]}    {'TXTL_UTR_UTR1_F'                }
%     {[  462.8714]}    {'TXTL_P70_RNAPbound_Kd'          }
%     {[    4.4817]}    {'TXTL_P70_RNAPbound_F'           }
%     {[   20.0374]}    {'TXTL_RNAPBOUND_TERMINATION_RATE'}
%     {[   19.0278]}    {'TXTL_NTP_RNAP_1_Kd'             }
%     {[         1]}    {'TXTL_NTP_RNAP_1_F'              }
%     {[1.1990e+06]}    {'TXTL_NTP_RNAP_2_Kd'             }
%     {[         1]}    {'TXTL_NTP_RNAP_2_F'              }
%     {[  703.8744]}    {'TL_AA_Kd'                       }
%     {[    0.7408]}    {'TL_AA_F'                        }
%     {[2.0007e+06]}    {'TL_AGTP_Kd'                     }
%     {[    0.3012]}    {'TL_AGTP_F'                      }
%     {[   38.2385]}    {'TXTL_RIBOBOUND_TERMINATION_RATE'}
%     {[6.5322e+06]}    {'TXTL_RNAdeg_Kd'                 }
%     {[         1]}    {'TXTL_RNAdeg_F'                  }
%     {[    0.1453]}    {'TXTL_RNAdeg_kc'                 }
%     {[    9.0194]}    {'RNAP'                           }
%     {[   44.3395]}    {'Ribo'                           }
%     {[3.2219e+04]}    {'RNase'                          }
    


% the parameters to set in the model are the parameters in the mcmc_info
% constructor file, then we can use the resulting model with the txtl_plot
% command to generate a prediction. 
%%
parvals = struct('paramNames', mi(2).namesOrd,...
    'paramVals', num2cell((marrayOrd(:, end, end))),...
    'reactionString', 'global'...
    );
% simulate model with that parameter point
Mobj = model_dsg2014_regen('plotmode', true, 'initialDNA', 30, 'paramInfo', parvals)
% plot trajectories for that model . 

doseStruct = struct('dosedNames', {{'DNA p70--utr1--deGFP'}}, ...
    'dosedVals', 30);

mcmc_plot_txtl(Mobj, parvals,doseStruct, ...
    'model_info', mi(2), 'master_info', mai) % there is only a single geometry, so the
% 'which_geometry' argument does not need to be specified. 

