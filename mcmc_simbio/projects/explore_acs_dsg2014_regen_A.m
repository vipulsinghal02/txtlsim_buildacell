% Explore all the data collected on AWS (m5, c5 instances, usually 4x, 
% 12x and 24x large) and on the local machine (Mac Book Pro, mid 2015, 
% quad core 15 inch, academic individual license)

%% Sim set 1: 24 core machine (m5.12xlarge)
% Explore the simulations from the 24 core machine, where we used 3200
% walkers.
tstamptouse1 = '10pct_20190203_173126';
tstamptouse2 = 'pt1pct_20190203_173126';
projdir = [pwd '/mcmc_simbio/projects/proj_acs_dsg2014_regen_A'];
marray = mcmc_get_walkers({tstamptouse1, tstamptouse2}, {1:16, 1:13}, projdir);

msubsamp = marray(:, 1:100:end, :);
% 
% load([projdir...
%     '/simdata_' tstamptouse1 '/full_variable_set_' tstamptouse1], 'mi',...
%     'mcmc_info', 'data_info', 'mai', 'ri')
%
% plot the trace plots
parnames = ...
    [{'TX_{cat}'     }
    {'\tau_{atp}'   }
    {'pol_{Kd}'     }
    {'pol_{kf}'     }
    {'pol_{term}'   }
    {'n_{Kd1}'      }
    {'n_{kf1}'      }
    {'n_{Kd2}'      }
    {'n_{kf2}'      }
    {'RNAse_{Kd}'   }
    {'RNAse_{kf}'   }
    {'RNAse_{cat}'  }
    {'pol'          }
    {'RNase'        }
    {'TL_{cat}'     }
    {'Ribo_{Kd}'    }
    {'Ribo_{kf}'    }
    {'aa_{Kd}'      }
    {'aa_{kf}'      }
    {'TL_n_{Kd}'    }
    {'TL_n_{kf}'    }
    {'Ribo_{term}'  }
    {'Ribo'         }];

%%
mcmc_plot(msubsamp(1:10, :,:), parnames(1:10),...
    'savematlabfig', false, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', ['20190203_173126']);

%% 
mfinal = marray(:,1:10:end,(end-20):end);
mfinal = mfinal(:,:)';
%%
mcmc_plot(mfinal(:,1:10), parnames(1:10),...
    'savematlabfig', false, 'savejpeg', false,...
    'projdir', projdir, 'tstamp', ['20190203_173126']);

% plot the autocorrelation plot


% plot the corner plot

% plot the trajectory fits.














%% Sim set 2: 48 core machine (m5.24xlarge)

tstamptouse = tstamp_appended;
marray = mcmc_get_walkers({tstamptouse}, {1:ri.nIter}, projdir);
mcmc_plot(marray([1 2 4], :,:), mai.estNames([1 2 4]),...
    'savematlabfig', false, 'savejpeg', true,...
    'projdir', projdir, 'tstamp', tstamptouse);

titls = {};
lgds = {};
mvarray = masterVecArray(marray, mai);
marrayOrd = mvarray(mi(1).paramMaps(mi(1).orderingIx),:,:);
mcmc_trajectories(mi(1).emo, di(mi(2).dataToMapTo), mi(1), marrayOrd,...
    titls, lgds,...
    'SimMode', 'meanstd', 'savematlabfig', false, 'savejpeg', true,...
    'projdir', projdir, 'tstamp', tstamptouse);

marrayOrd = mvarray(mi(2).paramMaps(mi(2).orderingIx),:,:);
titls = {};
mcmc_trajectories(mi(2).emo, di(mi(2).dataToMapTo), mi(2), marrayOrd,...
    titls, lgds,...
    'SimMode', 'meanstd', 'savematlabfig', false, 'savejpeg', true,...
    'projdir', projdir, 'tstamp', tstamptouse);

