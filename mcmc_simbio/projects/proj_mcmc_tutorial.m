%% First tutorial file for the mcmc_simbio package
% proj_mcmc_tutorial.m - Basic toturial of the mcmc_simbio package
% demonstrating the estimation of parameters for a constitutive gene
% expression circuit modeled as a single enzymatic reaction. 

%% Initializing the toolbox
% If you have not already initialized the txtlsim and mcmc_simbio
% toolboxes, initialize them by running the txtl_init and mcmc_init
% commands in the command line. You need your working directory to be the
% main directory where the txtlsim toolbox is stored (i.e., the directory
% in which directories like core, components, mcmc_simbio etc are stored).

%% Run project_init 
% This creates a directory within the projects directory where 
% the results of the simulation will be stored. The name of the directory
% will be the same as the name of this file (proj_mcmc_tutorial, in this 
% case). it also creates a timestamped subdirectory within this directory, 
% where the actual results are stored. If the top level directory
% (proj_mcmc_tutorial) already exists, then only the subdirectory is
% created. 
[tstamp, projdir, st] = project_init;

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
% mcmc_info_constgfp3i.m file to generate the mcmc_info struct that we need
% to define our parameter inference problem. 
% Please enter 'help mcmc_info' into the command window prompt to read more
% about this struct. Also, open the mcmc_info_constgfp3i file (edit
% mcmc_info_constgfp3i) to view how it is set up. 

mcmc_info = mcmc_info_constgfp3i(mobj);

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
cpol1 = 100; % nM

% Arrange the parameters in a log transformed vector. 
masterVector = log([rkfdG 
                    rkrdG
                    rkcp1
                    cpol1]);

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

di = data_artificial_v2({mobj}, {0:180:7200}, {mi.measuredSpecies}, ...
    {mi.dosedNames}, {mi.dosedVals}, {mi.namesUnord},...
    {exp(masterVector), [exp(masterVector(1:end-2)); 0.024; 200]});

da_extract1 = di(1).dataArray;
tv = di(1).timeVector;

%% Plot the artificial data
% we can plot the data using the mcmc_trajectories function. See its help
% file for usage information. 
mcmc_trajectories([], di, [], [], [], [], 'just_data_info', true);

%% Run the MCMC 
ri = mcmc_info.runsim_info;

mai = mcmc_info.master_info;

mi1 = mcmc_runsim_v2(tstamp, projdir, di(1), mcmc_info,...
   'InitialDistribution', 'LHS', 'multiplier', 2); % 'InitialDistribution', 'gaussian'

%%  plot stuff 
tstamptouse = tstamp; 
marray = mcmc_get_walkers({tstamptouse}, {1:ri.nIter}, projdir);
mcmc_plot(marray, mai.estNames, 'savematlabfig', true, 'savejpeg', true,...
    'projdir', projdir, 'tstamp', tstamptouse);
titls = {'dna 10'; 'dna 30';'dna 60'};
lgds = {};
mvarray1 = masterVecArray(marray, mai);
marrayOrd = mvarray1(mi1.paramMaps(mi1.orderingIx, 1),:,:);
fhandle = mcmc_trajectories(mi1.emo, di(1), mi1, marrayOrd, titls, lgds,...
    'SimMode', 'curves', 'savematlabfig', true, 'savejpeg', true,...
    'projdir', projdir, 'tstamp', tstamptouse);


%  Vipul Singhal, 
%  California Institute of Technology
%  2018
