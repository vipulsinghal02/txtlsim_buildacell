


% In this file we try to find a set of parameter values that allow the
% model from model_dsg2014_regen to fit the data from data_dsg2014_full.
% Vipul Singhal, 2019

% simulate txtl model with custom parameter values, and look at the species
% plots as specified by mcmc_info object.

% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_acs_dsg2014_regen_G'];
addpath(projdir)

jpgsave = true;
figsave = true;

% Load model, mcmc_info, and data_info.
mobj = model_dsg2014_regen;
mcmc_info = mcmc_info_dsg2014_regen_G(mobj);
di = data_dsg2014_full;

mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

% plot data from existing simulations.
tsIDtouse = 1;
plotflag = true;
switch tsIDtouse
    case 1 
        ts1 = '20190209_191801_1_127';
        ts2 = '20190210_030032_1_127';
        ts3 = '20190210_123436_1_32';
        ts4 = '20190210_185320_1_3';
        nIterID = {1:20 1:29 1:28 1:80};
        tstamp = {ts1 ts2 ts3};
        load([projdir '/simdata_' ts1temp '/full_variable_set_' ts1temp '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
    
end

mai.masterVector


marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray = marray_full(:, :, 1:end);
parnames = ...
    [{'TX_{cat}'    }
    {'pol_{Kd}'     }
    {'pol'          }
    {'TL_{cat}'     }
    {'Ribo_{Kd}'    }
     {'Ribo'         }];
%

if plotflag
    close all
    % Plot trace and corner (posterior distribution) plots
mcmc_plot(marray(:, 1:10:end,1:50:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', ts1, 'extrafignamestring', 'WithTransient');
mcmc_plot(marray(:, 1:3:end,(120*50):50:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', ts1, 'extrafignamestring', 'BurnedIn');
%     mcmc_plot(marray(6+(1:8), 1:20:end,1:end), parnames(6+(1:8)));

    figure
    [C,lags,ESS]=eacorr(marray(:, :,(120*50):end));%10000:end
    plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    grid on
    xlabel('lags')
    ylabel('autocorrelation');
    text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',...
        ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    title('Markov Chain Auto Correlation')
    %
    titls = arrayfun(@num2str, mi(1).dosedVals, 'UniformOutput', false);
    titls_array = cell(length(titls), 1, length(mi(1).measuredSpeciesIndex));
    ms = {'MG aptamer', 'deGFP'};
    for i = 1:length(mi(1).measuredSpeciesIndex)
        for j = 1:length(titls)
            titls_array(j, 1, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, Exp data']};
            titls_array(j, 2, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, MCMC samples']};
        end
    end
    mvarray = masterVecArray(marray, mai);
    samplePoints = ceil(size(mvarray, 3) * [.9, 1]);
    marrayOrd = mvarray(mi(1).paramMaps(mi(1).orderingIx),:,samplePoints);
    mcmc_trajectories(mi(1).emo, data_info(mi(1).dataToMapTo), mi(1), marrayOrd,...
        titls_array, {},...
        'SimMode', 'meanstd', 'separateExpSim', true);
end
marrayOrd(:,1:5,end)
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]
%