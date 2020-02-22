
% In this file we try to find a set of parameter values that allow the
% model from model_dsg2014_regen to fit the data from data_dsg2014_full.
% Vipul Singhal, 2019

% simulate txtl model with custom parameter values, and look at the species
% plots as specified by mcmc_info object.

% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_acs_dsg2014_regen_H'];
addpath(projdir)

jpgsave = true;
figsave = false;

% Load model, mcmc_info, and data_info.
mobj = model_dsg2014_regen;
mcmc_info = mcmc_info_dsg2014_regen_H(mobj);
di = data_dsg2014_full;

mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

% plot data from existing simulations.
tsIDtouse = 3;
plotflag = true;
switch tsIDtouse
    case 1
        ts1 = '20190211_122853_1_290';
        ts2 = '20190211_213614_1_29';
        ts3 = '20190213_141754_1_116';
        tstamp = {ts1 ts2 ts3};
        nIterID = {1:10 1:10 1:20};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
        
    case 2
        ts1 = '20190211_122853_1_290'; %0.005
        ts2 = '20190213_173314_1_174'; % 0.003
        ts3 = '20190213_173314_2_58'; % 0.001
        ts4 = '20190213_173314_3_52'; %0.009
        ts5 = '20190213_173314_4_41'; % 0.007
        ts6 = '20190213_173314_5_29'; % 0.005
        
    tstamp = {ts1 ts2 ts3 ts4 ts5 ts6};
        nIterID = {1:10 1:5 1:5 1:5 1:5 1:5};
        load([projdir '/simdata_' ts6 '/full_variable_set_' ts6 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
    case 3
        ts1 = '20190211_122853_1_290'; %0.005
        ts2 = '20190215_102557_1_232';
        ts3 = '20190215_102557_2_174';
        ts4 = '20190215_102557_3_116';
        
    tstamp = {ts1 ts2 ts3 ts4};
        nIterID = {1:10 1:20 1:20 1:20};
        load([projdir '/simdata_' ts2 '/full_variable_set_' ts2 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
end

mai.masterVector


marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray = marray_full(:,:,1:end);
parnames = ...
    [{'TX_{cat}'    }
    {'\tau_{atp}'   }
    {'\delta_{atp}' }
    {'\alpha_{atp}' }
    {'pol_{Kd}'     }
    {'pol_{term}'   }
    {'n_{Kd1}'      }
    {'n_{Kd2}'      }
    {'RNAse_{Kd}'   }
    {'RNAse_{cat}'  }
    {'pol'          }
    {'RNase'        }
    {'TL_{cat}'     }
    {'GFP_{mat}'    }
    {'Ribo_{Kd}'    }
    {'aa_{Kd}'      }
    {'TL_n_{Kd}'    }
    {'Ribo_{term}'  }
    {'Ribo'         }];
% 

% 
% 
%     {'TX_{cat}'    }    {'TX_elong_glob'                  }
%     {'\tau_{atp}'  }    {'AGTPdeg_time'                   }
%     {'\delta_{atp}'}    {'AGTPdeg_rate'                   }
%     {'\alpha_{atp}'}    {'AGTPreg_ON'                     }
%     {'pol_{Kd}'    }    {'TXTL_P70_RNAPbound_Kd'          }
%     {'pol_{term}'  }    {'TXTL_RNAPBOUND_TERMINATION_RATE'}
%     {'n_{Kd1}'     }    {'TXTL_NTP_RNAP_1_Kd'             }
%     {'n_{Kd2}'     }    {'TXTL_NTP_RNAP_2_Kd'             }
%     {'RNAse_{Kd}'  }    {'TXTL_RNAdeg_Kd'                 }
%     {'RNAse_{cat}' }    {'TXTL_RNAdeg_kc'                 }
%     {'pol'         }    {'RNAP'                           }
%     {'RNase'       }    {'RNase'                          }
%     {'TL_{cat}'    }    {'TL_elong_glob'                  }
%     {'GFP_{mat}'   }    {'TXTL_PROT_deGFP_MATURATION'     }
%     {'Ribo_{Kd}'   }    {'TXTL_UTR_UTR1_Kd'               }
%     {'aa_{Kd}'     }    {'TL_AA_Kd'                       }
%     {'TL_n_{Kd}'   }    {'TL_AGTP_Kd'                     }
%     {'Ribo_{term}' }    {'TXTL_RIBOBOUND_TERMINATION_RATE'}
%     {'Ribo'        }    {'Ribo'                           }
%
%%
if plotflag
    close all
    % Plot trace and corner (posterior distribution) plots
%     mcmc_plot(marray(:, 1:3:end,1:50:end), parnames(:));
%     mcmc_plot(marray(:, 1:end,5000:5:end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', ts1, 'extrafignamestring', 'BurnedInAllWalkers');

% 
%     mcmc_plot(marray(:, 1:20:end,5000:10:end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', ts1, 'extrafignamestring', 'BurnedIn');

    mcmc_plot(marray(:, 1:end,(end-2000):10:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', ts2, 'extrafignamestring', 'last300');

    mcmc_plot(marray(:, 1:30:end,1:10:end), parnames(:),...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', ts2, 'extrafignamestring', 'WithTransient');
% 
%     mcmc_plot(marray(:, 1:end,(end-4000):50:end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', ts1, 'extrafignamestring', 'WithTransient');
% 
%     mcmc_plot(marray(:, 1:40:end,(end-4000):10:end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', jpgsave,...
%     'projdir', projdir, 'tstamp', ts1, 'extrafignamestring', 'WithTransient');
%     figure
%     [C,lags,ESS]=eacorr(marray(:, :,5000:end));%10000:end
%     plot(lags,C,'.-',lags([1 end]),[0 0],'k');
%     grid on
%     xlabel('lags')
%     ylabel('autocorrelation');
%     text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',...
%         ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
%     title('Markov Chain Auto Correlation')
%     
%%    
    mvarray = masterVecArray(marray, mai);
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
    mcmc_trajectories(mi(1).emo, data_info(mi(1).dataToMapTo), mi(1), marrayOrd,...
        titls_array, {},...
        'SimMode', 'meanstd', 'separateExpSim', true,...
        'savematlabfig', figsave, 'savejpeg', jpgsave,...
        'projdir', projdir, 'tstamp', ts2,...
        'extrafignamestring', 'RNAspike');
    
    marrayOrd = mvarray(mi(2).paramMaps(mi(2).orderingIx),:,samplePoints);
    titls = arrayfun(@num2str, mi(2).dosedVals, 'UniformOutput', false);
    titls_array = cell(length(titls), 2, length(mi(2).measuredSpeciesIndex));
    ms = {'MG aptamer', 'deGFP'};
    for i = 1:length(mi(2).measuredSpeciesIndex)
        for j = 1:length(titls)
            titls_array(j, 1, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, Exp data']};
            titls_array(j, 2, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, MCMC samples']};
        end
    end
    mcmc_trajectories(mi(2).emo, data_info(mi(2).dataToMapTo), mi(2), marrayOrd,...
        titls_array, {},...
        'SimMode', 'meanstd', 'separateExpSim', true,...
        'savematlabfig', figsave, 'savejpeg', jpgsave,...
        'projdir', projdir, 'tstamp', ts2,...
        'extrafignamestring', 'MGa_deGFP');
    %%
end
marrayOrd(:,1:5,end)
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]
