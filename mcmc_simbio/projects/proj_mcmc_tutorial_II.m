%% Second tutorial for the mcmc_simbio package
% proj_mcmc_tutorial_II.m - tutorial of the mcmc_simbio package
% demonstrating the estimation of parameters for a constitutive gene
% expression circuit modeled as a single enzymatic reaction. In this 
% example, we demonstrate the concurrent parameter
% estimation capability, where we have the constitutive expression circuit
% in two environments, with each environment having its own environment
% specific parameters (ESPs), while the circuit having a single set of
% circuit specific parameters. 

%% Initializing the toolbox
% If you have not already initialized the txtlsim and mcmc_simbio
% toolboxes, initialize them by running the txtl_init and mcmc_init
% commands in the command line. You need your working directory to be the
% main directory where the txtlsim toolbox is stored (i.e., the directory
% in which directories like core, components, mcmc_simbio etc are stored).

%% Run project_init 
% This creates a directory within the projects directory where 
% the results of the simulation will be stored. The name of the directory
% will be the same as the name of this file (proj_mcmc_tutorial_II, in this 
% case). it also creates a timestamped subdirectory within this directory, 
% where the actual results are stored. If the top level directory
% (proj_mcmc_tutorial_II) already exists, then only the subdirectory is
% created. 
delete(gcp('nocreate'))
parpool(47)
[tstamp, projdir, st] = project_init;
prevtstamp = '20190131_181526'
%% Define the MATLAB Simbiology model 
% We use the file model_protein3.m to define a constitutive gene expression 
% model using a single enzymatic step. The reactions and species that it
% sets up are: 
% 
% dG + pol <-> dG_pol  (k_f, k_r) 
% dG_pol -> dG + pol + pG (kc) 

mobj = model_protein3;

% The species of the model can be visualized as follows: 
mobj.species

% The reactions may be visualized as
mobj.reactions
% For more on MATLAB Simbiology, see the Simbiology 
% <https://www.mathworks.com/help/simbio/gs/simbiology-command-line-tutorial.html
% reference> page. 

%% Defining the experiment / model arrangement. 
% We can define the experimental setup and how it related to data, the
% Simbiology model and the estimation problem using what we call an
% mcmc_info struct. For this example, we will be using an
% mcmc_info_constgfp3ii.m file to generate the mcmc_info struct that we need
% to define our parameter inference problem. 
% Please enter 'help mcmc_info' into the command window prompt to read more
% about this struct. Also, open the mcmc_info_constgfp3ii file (enter 'edit
% mcmc_info_constgfp3ii' into the command prompt) to learn how it is set up. 

mcmc_info = mcmc_info_constgfp3ii(mobj);

%% Creating artificial data to fit the model to. 
% Instead of using real data, we will create artificial data for
% demonstration purposes. We will use the data_artificial_v2 fucntion to do
% this. 

% Get the model_info struct needed to generate the artificial data 
mi = mcmc_info.model_info; 

% A list of nominal parameter values to use to generate the data. 
rkfdG = 5; % nM-1s-1
rkrdG = 300; % s-1
rkcp1 = 0.012; %s-1
rkcp2 = 0.024; %s-1
cpol1 = 100; % nM
cpol2 = 200; % nM

% Arrange the parameters in a log transformed vector. 
masterVector = log([...
rkfdG 
rkrdG
rkcp1
rkcp2
cpol1
cpol2]);

% Supply the experimental setup information to the data_artificial_v2 
% function so that it can generate the data_info struct that contains the 
% artificial data. 
% type 'help data_artificial_v2' into the command window prompt to read
% more about this function. For our purposes we simply note that we need to
% specify our Simbiology model object, a set of timepoints to report the
% output trajectories for, the list of measured species' names for our 
% model, the list of dosed species' names, the matrix of dosed values, the
% names of the species and parameters to set values for in the model
% (namesUnord), and the non-log-transformed values as a vector. All of
% these arguments must be encapsulated in cells. 

di = data_artificial_v2({mobj}, {0:180:7200}, {mi.measuredSpecies},...
    {mi.dosedNames}, {mi.dosedVals}, {mi.namesUnord},...
     {exp(masterVector([1:2 3 5])), exp(masterVector([1:2 4 6]))});

 
da_extract1 = di(1).dataArray;
da_extract2 = di(2).dataArray;
tv = di(1).timeVector;

%% Plot the artificial data
% we can plot the data using the mcmc_trajectories function. See its help
% file for usage information. 
%mcmc_trajectories([], di, [], [], [], [], 'just_data_info', true);

%% Run the MCMC 
ri = mcmc_info.runsim_info;

mai = mcmc_info.master_info;


specificprojdir = [projdir '/simdata_' prevtstamp];

% load mcmc_info    and the updated model_info
SS = load([specificprojdir '/full_variable_set_' prevtstamp], 'mcmc_info');

marray = mcmc_get_walkers({prevtstamp}, {SS.mcmc_info.runsim_info.nIter},...
    projdir);
% assume the projdir where this data is stored is the same one as the
% one created at the start of this file

%%
pID = 1:length(mai.estNames);
marray_cut = mcmc_cut(marray, pID, flipud((mai.paramRanges)'));
if size(marray_cut, 2) < ri.nW
    warning('too few initial points, using a few timesteps from previous runs to create initial walker positions.');
    walker_timepoints = ceil(linspace(ceil(size(marray_cut,3))/2, size(marray_cut,3), ceil(ri.nW/size(marray_cut, 2))))
    minit = marray_cut(:,:, walker_timepoints(1));
    for i = 2:length(walker_timepoints)
        minit = [minit marray_cut(:,:,walker_timepoints(i)) ];
    end
    minit = minit(:, 1:ri.nW);
else % there are enough points, just pick the number needed. 
    minit = marray_cut(:,1:ri.nW,end);
end

%%

% now run the simulation. 

mi = mcmc_runsim_v2(tstamp, projdir, di, mcmc_info,...
   'UserInitialize', minit, 'multiplier', 2,...
   'pausemode', false); 
% 
% mi = mcmc_runsim_v2(tstamp, projdir, di, mcmc_info,...
%    'InitialDistribution', 'LHS', 'multiplier', 2,...
%    'pausemode', false); 
% 'InitialDistribution', 'gaussian'
%  'UserInitialize', marray_cut(:,:,end)

%%  plot stuff 



% load([pwd '/mcmc_simbio/projects/proj_mcmc_tutorial_II/'...
%     'simdata_20180802_221756/full_variable_set_20180802_221756'])
%%
%  projdir = [pwd '/mcmc_simbio/projects/proj_mcmc_tutorial_II']
% di = data_info
%  tstamptouse = tstamp; 
% marray = mcmc_get_walkers({tstamptouse}, {1:ri.nIter}, projdir);
% mcmc_plot(marray([1 2 4], :,:), mai.estNames([1 2 4]),...
%     'savematlabfig', true, 'savejpeg', true,...
%     'projdir', projdir, 'tstamp', tstamptouse,...
%     'extrafignamestring', '_extract1');
% %
% figure
% mcmc_plot(marray([1 3 5], :,:), mai.estNames([1 3 5]),...
%     'savematlabfig', true, 'savejpeg', true,...
%     'projdir', projdir, 'tstamp', tstamptouse,...
%     'extrafignamestring', '_extract2');
% titls = {'E1 dG 10';'E1 dG 30';'E1 dG 60';};
% lgds = {};
% mvarray = masterVecArray(marray, mai);
% marrayOrd = mvarray(mi(1).paramMaps(mi(1).orderingIx, 1),:,:);
% fhandle = mcmc_trajectories(mi(1).emo, di(1), mi(1), marrayOrd,...
%     titls, lgds,...
%     'SimMode', 'meanstd', 'savematlabfig', true, 'savejpeg', true,...
%     'projdir', projdir, 'tstamp', tstamptouse, 'extrafignamestring',...
%     '_extract1');
% marrayOrd = mvarray(mi(1).paramMaps(mi(1).orderingIx, 2),:,:);
% titls = {'E2 dG 10';'E2 dG 30';'E2 dG 60';};
% fhandle = mcmc_trajectories(mi(1).emo, di(2), mi(1), marrayOrd,...
%     titls, lgds,...
%     'SimMode', 'meanstd', 'savematlabfig', true, 'savejpeg', true,...
%     'projdir', projdir, 'tstamp', tstamptouse, 'extrafignamestring',...
%     '_extract2');

% 
% load([pwd '/mcmc_simbio/projects/proj_mcmc_tutorial_II/'...
%     'simdata_20190131_181526/full_variable_set_20190131_181526'])
%%


%  Vipul Singhal, 
%  California Institute of Technology
%  2018
