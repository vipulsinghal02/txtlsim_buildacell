% Plotting the core mechanism characterization results for the paper:
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
% You need around 8 GB
% of free hard disk space to download the full dataset. Download at: 
% https://www.dropbox.com/s/alx379b65ga40mh/core_characterization_files.zip?dl=0
% Extract the project folders and place them in the projects directory. 
% See lines 55, 85, 117, etc for examples of what the final path should
% look like. 
% 3. Make sure you have 8 GB of free ram. Otherwise, just run each case 
% (the for loop around line 41 below) separately, and the RAM usage should 
% never go over 3GB or so. 
% 4. The plots should be saved in the a subdirectory of the respective 
% project directory. The subdirectory is given by the final timestamp of that
% case. So for example, for the runN = 0 case, 
% the project directory is 
% ./mcmc_simbio/projects/proj_vnprl (line 55)
% and the subdirectory is 
% simdata_20200513_063545_4_19 (line 66)

close all 
clear all
figsave = false;
jpgsave = true;
txtldir = txtl_init;
mcmc_init;
% ----------------------------------------------------------------------- %
% Select Run number 
% Change [0 1 2 3 4 5] to [0] to run a single run (run 0 in this case) at a time
% This way, instead of 8GB of RAM, you only need 3 GB of RAM.
% ----------------------------------------------------------------------- %
for runN = [0 1 2 3 4 5] % or 0 1 2 3 4 5
    % note that run 0 is a special run that has the same parameters as run 3,
    % but uses composite data.
    % ----------------------------------------------------------------------- %
    % Load the models and the data
    % ----------------------------------------------------------------------- %
    % H = run 1
    % I = run 2
    % F = run 3
    % D = run 4
    % G = run 5
    % set the run specific information.
    switch runN
        case 0
            projdir = './mcmc_simbio/projects/proj_vnprl';
            mobj = model_dsg2014_regen;
            mcmc_info = mcmc_info_vnprl_F3(mobj);
            di = data_VNPRL2011;
            ts14 = '20200510_045912_1_95'; %300
            ts15 = '20200510_045912_2_48'; %600
            ts16 = '20200512_192355_1_43';
            ts17 = '20200512_192355_2_38';
            ts18 = '20200513_063545_1_33';
            ts19 = '20200513_063545_2_29';
            ts20 = '20200513_063545_3_24';
            ts21 = '20200513_063545_4_19';
            tstamp = {ts14 ts15 ts16 ts17 ts18 ts19 ts20 ts21}; %
            nIterID = {1:30 1:30 1:10  1:10 1:15 1:15 1:15 1:15};%
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
            parorder= [2 3 1 4 5 8 10 11 12 13 7 6 9];
        case 1 % alt code: H
            projdir = './mcmc_simbio/projects/proj_acs_dsg2014_regen_H';
            mobj = model_dsg2014_regen;
            mcmc_info = mcmc_info_dsg2014_regen_H(mobj);
            di = data_dsg2014_full;
            ts1 = '20190211_122853_1_290';
            ts2 = '20190215_102557_1_232';
            ts3 = '20190215_102557_2_174';
            ts4 = '20190215_102557_3_116';
            tstamp = {ts1 ts2 ts3 ts4}; %
            nIterID = {1:10 1:20 1:20 1:20};%
            parnames = ...
                [{'TX_{cat}'    }%1
                {'\tau_{atp}'   }
                {'\delta_{atp}' }
                {'\alpha_{atp}' }%4
                {'pol_{Kd}'     }
                {'pol_{term}'   }%6
                {'n_{Kd1}'      }
                {'n_{Kd2}'      }%8
                {'RNAse_{Kd}'   }
                {'RNAse_{cat}'  }
                {'pol'          }%11
                {'RNase'        }
                {'TL_{cat}'     }% 13
                {'GFP_{mat}'    }%14
                {'Ribo_{Kd}'    }
                {'aa_{Kd}'      }
                {'TL_n_{Kd}'    }%17
                {'Ribo_{term}'  }
                {'Ribo'         }];%19
            parorder = [2 3 4 1 5 6 11 7 8 13 15 18 19 16 17 10 9 12 14 ];
        case 2 %alt code: I
            projdir = './mcmc_simbio/projects/proj_acs_dsg2014_regen_I';
            mobj = model_dsg2014_regen;
            mcmc_info = mcmc_info_dsg2014_regen_I(mobj);
            di = data_dsg2014_full;
            ts1 = '20190211_123400_1_290';
            ts2 = '20190211_214449_1_29';
            tstamp = {ts1 ts2}; %
            nIterID = {1:9 1:20};%
            parnames = ...
                [{'TX_{cat}'    }% 1 
                {'\tau_{atp}'   }
                {'\delta_{atp}' }
                {'\alpha_{atp}' } %4
                {'pol_{Kd}'     }
                {'pol_{term}'   }
                {'RNAse_{Kd}'   }
                {'RNAse_{cat}'  } %8 
                {'pol'          }
                {'RNase'        }
                {'TL_{cat}'     }
                {'GFP_{mat}'    }%12
                {'Ribo_{Kd}'    }
                {'Ribo_{term}'  }
                {'Ribo'         }];
            
            parorder = [2 3 4 1 5 6 9 11 13 14 15 8 7 10 12];
        case 3 % alt code: F1
            projdir = './mcmc_simbio/projects/proj_acs_dsg2014_regen_F1';
            mobj = model_dsg2014_regen;
            mcmc_info = mcmc_info_dsg2014_regen_F1(mobj);
            di = data_dsg2014_full;
            ts1 = '20190210_184039_1_195';
            ts2 = '20190210_184039_2_20';
            ts3 = '20190211_054328_1_20';
            ts4 = '20190211_140024_1_20';
            ts5 = '20190212_053457_1_20';
            ts6 = '20190212_114314_1_20';
            tstamp = {ts6};%ts1 ts2 ts3 ts4 ts5 
            nIterID = {1:22};%1:10 1:10 1:15 1:27 1:11 
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
            parorder= [2 3 1 4 5 8 10 11 12 13 7 6 9];
        case 4% D
            projdir = './mcmc_simbio/projects/proj_acs_dsg2014_regen_D';
            mobj = model_dsg2014_regen;
            mcmc_info = mcmc_info_dsg2014_regen_D(mobj);
            di = data_dsg2014_full;
            ts1 = '20190209_184946_1_44';
            ts2 = '20190210_072127_1_20';
            tstamp = {ts1 ts2};
            nIterID = {1:24 1:16};
            parnames = ...
                [{'TX_{cat}'    }
                {'pol_{Kd}'     }
                {'pol_{term}'   }
                {'RNAse_{cat}'  }
                {'pol'          }
                {'RNase'        }
                {'TL_{cat}'     }
                {'Ribo_{Kd}'    }
                {'Ribo_{term}'  }
                {'Ribo'         }];
            
            parorder= [1 2 3 5 7 8 9 10 4 6];
        case 5
            projdir = './mcmc_simbio/projects/proj_acs_dsg2014_regen_G';
            mobj = model_dsg2014_regen;
            mcmc_info = mcmc_info_dsg2014_regen_G(mobj);
            di = data_dsg2014_full;
            ts1 = '20190209_191801_1_127';
            ts2 = '20190210_030032_1_127';
            ts3 = '20190210_123436_1_32';
            ts4 = '20190210_185320_1_3';
            nIterID = {1:20 1:40 1:28};% 1:80
            tstamp = {ts1 ts2 ts3};% ts4
            parnames = ...
                [{'TX_{cat}'    }
                {'pol_{Kd}'     }
                {'pol'          }
                {'TL_{cat}'     }
                {'Ribo_{Kd}'    }
                {'Ribo'         }];
            
            
            parorder= [1 2 3 4 5 6];
    end
    
    
    addpath(projdir)
    mi = mcmc_info.model_info;
    ri = mcmc_info.runsim_info;
    mai = mcmc_info.master_info;
    load([projdir '/simdata_' tstamp{end} '/full_variable_set_' tstamp{end} '.mat'], ...
        'mi',...
        'mcmc_info', 'data_info', 'mai', 'ri');
    tsToSave = tstamp{end};
    mai.masterVector
    marray = mcmc_get_walkers(tstamp,nIterID, projdir);
    
    %% Plot the results
    
    switch runN
        case 0
            mcmc_plot(marray(parorder, 1:5:end,1:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', true,...
                'plotDistribution', false, ...
                'fontsize', 16,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['all_iter_run' num2str(runN)]);
            
            % last 4k iterations are burned in.
            mcmc_plot(marray(parorder, 1:end,(end-200):50:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', false,'fontsize', 16,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['burnedin_run' num2str(runN)]);
        case 1
            % note that the original data is not thinned. Usually we thin the
            % iterations by 10. Here, due to the sheer size of the data, we
            % thin by 50
            mcmc_plot(marray(parorder, 1:end,1:70:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', true,...
                'plotDistribution', false, ...
                'fontsize', 14,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['all_iter_run' num2str(runN)]);
            
            mcmc_plot(marray(parorder, 1:end,(end-4000):1000:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', false,'fontsize', 14,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['burnedin_run' num2str(runN)]);
        case 2
            mcmc_plot(marray(parorder, 1:5:end,1:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', true,...
                'plotDistribution', false, ...
                'fontsize', 14,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['all_iter_run' num2str(runN)]);
            
            mcmc_plot(marray(parorder, 1:end,(end-1500):500:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', false,'fontsize', 14,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['burnedin_run' num2str(runN)]);
        case 3
            mcmc_plot(marray(parorder, 1:5:end,1:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', true,...
                'plotDistribution', false, ...
                'fontsize', 16,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['all_iter_run' num2str(runN)]);
            
            mcmc_plot(marray(parorder, 1:end,(end-3000):500:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', false,'fontsize', 16,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['burnedin_run' num2str(runN)]);
        case 4
            mcmc_plot(marray(parorder, 1:2:end,1:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', true,...
                'plotDistribution', false, ...
                'fontsize', 16,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['all_iter_run' num2str(runN)]);
            
            % last 4k iterations are burned in.
            mcmc_plot(marray(parorder, 1:end,(end-1000):500:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', false,'fontsize', 16,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['burnedin_run' num2str(runN)]);
        case 5
            aa = 1;
            mcmc_plot(marray(parorder, 1:3:end,1:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', true,...
                'plotDistribution', false, ...
                'fontsize', 20,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['all_iter_run' num2str(runN)]);
            
            mcmc_plot(marray(parorder, 1:end,(end-2500):250:end), parnames(parorder),...
                'savematlabfig', figsave, 'savejpeg', jpgsave,...
                'plotChains', false,'fontsize', 20,...
                'projdir', projdir, 'tstamp', tsToSave, ...
                'extrafignamestring', ['burnedin_run' num2str(runN)]);
    end
    
    %%
    % ----------------------------------------------------------------------- %
    % Plot the trajectories.
    % ----------------------------------------------------------------------- %
    if runN ~= 5
        workingDir = [pwd '/mcmc_simbio/exp_data/public_data/'];
        % run merger
        %all_data_merger;
        load([workingDir 'mergedExperimentFiles.mat'])
        colorCodes = {'r','b','g','c','m','k',[1 1 .5],[.7 .5 .2],[0 1 .2],[.35 .8 .8],[.9 0 .4],[1 .2 .2]};
        
        RFUumConvert = 1.723;   % (1.723 a.u. = 1 nM)
        RFUumConvertMG = 7.75;  % (7.75 a.u. = 1 nM)
        if runN == 0
            ylims = [1000 50 12000];
        else
            ylims = [1000 50 12000].*[1 10 1.8];
        end
        lengthToPlotArray = [21, 76, 76];
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
        %
        mvarray = masterVecArray(marray(:,:,(end-300):100:end), mai);
        samplePoints = ceil(size(mvarray, 3) * [.9, 1]);
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
            [da{miID}, idxnotused{miID}] = ...
                simulatecurves(currmi.emo,marrayOrd(:,:)', 50, dose, tv, currmi.measuredSpecies);
        end
        
        t_vec_mins=mergedExpFile(1).t_vec/60;
        
        figure
        set(0,'Units','normalized')
        set(gcf,'Units', 'normalized')
        set(gcf, 'Position', [0.05, 0.1, 0.4, 0.6])
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
                
                if runN == 0
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
                else
                    for i=1:11
                        mgNoBgMean(:,i)=mean([ ...
                            expFile(1).Data(:,12+i,1)-mgBackground(:,i) ...
                            expFile(2).Data(:,12+i,1)-mgBackground(:,i) ...
                            expFile(2).Data(:,35+i,1)-mgBackground(:,i) ...
                            ],2)/RFUumConvertMG;
                        
                        mgNoBgErr(:,i)=std([ ...
                            expFile(1).Data(:,12+i,1)-mgBackground(:,i) ...
                            expFile(2).Data(:,12+i,1)-mgBackground(:,i) ...
                            expFile(2).Data(:,35+i,1)-mgBackground(:,i) ...
                            ],0,2)/sqrt(3)/RFUumConvertMG;
                    end
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
                %h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,plotSelect+11,1)
                %         /RFUumConvertMG,mergedExpFile(1).noBg_std(:,plotSelect+11)/RFUumConvertMG,0.3,colorCodes);
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
                
                
                if runN == 0
                    h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,11+plotSelect,2)/RFUumConvert/1.8,...
                        mergedExpFile(1).noBg_std(:,11+plotSelect,2)/RFUumConvert/sqrt(3)/1.8,0.3,colorzCell);
                else
                    h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,11+plotSelect,2)/RFUumConvert,...
                        mergedExpFile(1).noBg_std(:,11+plotSelect,2)/RFUumConvert/sqrt(3),0.3,colorzCell);
                end
                
                axis([0 (lengthToPlot-1)*timeinterval(count),0,ylims(count)])
                title([titleArray{count} ' (Exp)'])
                xlabel('time, minutes')
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
                
                [linehandle(dID), ptchhandle] =...
                    boundedline(tv(1:lengthToPlot)/60,...
                    mean(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID), 3),...
                    std(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID)...
                    +eps*randn(lengthToPlot, 1, 50), 0,3));
                set(ptchhandle, ...
                    'FaceColor', colorz(dID+2, :),...
                    'FaceAlpha', 0.25);
                set(linehandle(dID), ...
                    'Color', colorz(dID+2, :),...
                    'LineWidth', 1.5);
                
                hold on
            end
            grid on
            title([titleArray{count} ' (Fit)'])
            
            legend(linehandle(dosesToPlot{count}{:}), legends{count, :}, 'Location', legLoc{count}, 'FontSize', 14)
            legend('boxoff')
            axis([0 (lengthToPlot-1)*timeinterval(count) 0 ylims(count)])
            if count ==length(miToUse)
                xlabel('time, minutes', 'FontSize', 14)
            end
            ax = gca;
            ax.FontSize = 14;
        end
        specificproj = [projdir '/simdata_' tsToSave];
%         print(gcf, '-djpeg', '-r600', [specificproj '/traj' tsToSave])
%         print(gcf, '-dpng', '-r600', [specificproj '/traj' tsToSave])
        print(gcf, '-djpeg', '-r200', [specificproj '/traj_lowres_' tsToSave])
        print(gcf, '-dpng', '-r200', [specificproj '/traj_lowres_' tsToSave])        
        
    else
        
        workingDir = [pwd '/mcmc_simbio/exp_data/public_data/'];
        % run merger
        %all_data_merger;
        load([workingDir 'mergedExperimentFiles.mat'])
        
        colorCodes = {'r','b','g','c','m','k',[1 1 .5],[.7 .5 .2],[0 1 .2],[.35 .8 .8],[.9 0 .4],[1 .2 .2]};
        
        RFUumConvert = 1.723;   % (1.723 a.u. = 1 nM)
        RFUumConvertMG = 7.75;  % (7.75 a.u. = 1 nM)
        
        ylims = [50 12000].*[10 1.8];
        
        lengthToPlotArray = [76, 76];
        miToUse = [1 1]; % mi 1 has one plot (RNA) and mi 2 has 2 plots (RNA and protein)
        msToPlot = [1 2];
        timeinterval = [8 8];
        %
        titleArray = {'mRNA expression';
            'deGFP expression'};
        legLoc = {'NorthEast','NorthWest'};
        
        legends = [{{'0.5nM', '2nM', '5nM', '20nM'}};
            {{'0.5nM', '2nM', '5nM', '20nM'}};];
        dosesToPlot = {{[1 2 3 4]}
            {[1 2 3 4]}};
        
        %     clear marray
        yLab = {'RNA, nM', 'deGFP, nM'}
        mvarray = masterVecArray(marray(:,:,(end-300):100:end), mai);
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
            
            [da{miID}, idxnotused{miID}] = ...
                simulatecurves(currmi.emo,marrayOrd(:,:)', 50, dose, tv, currmi.measuredSpecies);
        end
        
        t_vec_mins=mergedExpFile(1).t_vec/60;
        
        figure
        set(0,'Units','normalized')
        set(gcf,'Units', 'normalized')
        set(gcf, 'Position', [0.05, 0.1, 0.4, 0.6])
        for count = 1:length(miToUse)  %1:length(mi)%1:
            colorz = flipud(parula(max(dosesToPlot{count}{:})+2));
            currmi = mi(miToUse(count));
            currdi = data_info(currmi.dataToMapTo);
            tv = currdi.timeVector;
            lengthToPlot = lengthToPlotArray(count);
            % experimental data
            subplot(2, 2, (count-1)*2+1)
            if count == 1 % for the MG aptamer and the deGFP, use the scripts provided by Zoltan.
                %%%%%  MG KINETICS WITH CORRECT BACKGROUNDS SUBTRACTED
                plotSelect=[ 5, 7, 8, 10];
                mgBackground = mergedExpFile(1).Data_mean(:,1:11,1);
                mgNoBgMean=mergedExpFile(1).Data_mean(:,12:22,1)-mgBackground;
                
                
                
                for i=1:11
                    mgNoBgMean(:,i)=mean([ ...
                        expFile(1).Data(:,12+i,1)-mgBackground(:,i) ...
                        expFile(2).Data(:,12+i,1)-mgBackground(:,i) ...
                        expFile(2).Data(:,35+i,1)-mgBackground(:,i) ...
                        ],2)/RFUumConvertMG;
                    
                    mgNoBgErr(:,i)=std([ ...
                        expFile(1).Data(:,12+i,1)-mgBackground(:,i) ...
                        expFile(2).Data(:,12+i,1)-mgBackground(:,i) ...
                        expFile(2).Data(:,35+i,1)-mgBackground(:,i) ...
                        ],0,2)/sqrt(3)/RFUumConvertMG;
                end
                
                
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
                %h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,plotSelect+11,1)
                %         /RFUumConvertMG,mergedExpFile(1).noBg_std(:,plotSelect+11)/RFUumConvertMG,0.3,colorCodes);
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
            elseif count == 2
                plotSelect=[5, 7, 8, 10];
                %         figure('Name','GFP kinetics');
                hold on
                for dID = dosesToPlot{count}{:}
                    colorzCell{dID} = colorz(dID+2, :);
                end
                
                
                if runN == 0
                    h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,11+plotSelect,2)/RFUumConvert/1.8,...
                        mergedExpFile(1).noBg_std(:,11+plotSelect,2)/RFUumConvert/sqrt(3)/1.8,0.3,colorzCell);
                else
                    h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,11+plotSelect,2)/RFUumConvert,...
                        mergedExpFile(1).noBg_std(:,11+plotSelect,2)/RFUumConvert/sqrt(3),0.3,colorzCell);
                end
                
                axis([0 (lengthToPlot-1)*timeinterval(count),0,ylims(count)])
                title([titleArray{count} ' (Exp)'])
                xlabel('time, minutes')
                ylabel('deGFP (uM)')
                grid on
                legendStr = cellstr(num2str(mergedExpFile(1).concentrations(plotSelect)', 'pr-gfp-mg15 %0.2g nM'));
                ax = gca;
                ax.FontSize = 14;
                %legend(h,legendStr,'Location','NorthWest')
            end
            
            % simulation results
            subplot(2, 2, (count-1)*2+2)
            
            for dID = dosesToPlot{count}{:}
                
                [linehandle(dID), ptchhandle] =...
                    boundedline(tv(1:lengthToPlot)/60,...
                    mean(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID), 3),...
                    std(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID)...
                    +eps*randn(lengthToPlot, 1, 50), 0,3));
                set(ptchhandle, ...
                    'FaceColor', colorz(dID+2, :),...
                    'FaceAlpha', 0.25);
                set(linehandle(dID), ...
                    'Color', colorz(dID+2, :),...
                    'LineWidth', 1.5);
                
                %
                %         plot(tv(1:lengthToPlot)/60,...
                %             mean(da{miToUse(count)}(1:lengthToPlot, msToPlot(count), :, dID), 3),...
                %             'LineWidth', 1.5,...
                %             'Color', colorz(dID+2, :))
                hold on
            end
            grid on
            title([titleArray{count} ' (Fit)'])
            
            legend(linehandle(dosesToPlot{count}{:}), legends{count, :},...
                'Location', legLoc{count}, 'FontSize', 14)
            legend('boxoff')
            axis([0 (lengthToPlot-1)*timeinterval(count) 0 ylims(count)])
            if count ==length(miToUse)
                xlabel('time, minutes', 'FontSize', 14)
            end
            ax = gca;
            ax.FontSize = 14;
        end
        specificproj = [projdir '/simdata_' tsToSave];
%         print(gcf, '-djpeg', '-r600', [specificproj '/traj' tsToSave])
%         print(gcf, '-dpng', '-r600', [specificproj '/traj' tsToSave])
        print(gcf, '-djpeg', '-r200', [specificproj '/traj_lowres_' tsToSave])
        print(gcf, '-dpng', '-r200', [specificproj '/traj_lowres_' tsToSave])
%         print(gcf, '-depsc', [specificproj '/traj_eps_' tsToSave])    
        
        
    end
end
