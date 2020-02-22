%% Second tutorial for the mcmc_simbio package
% proj_mcmc_tutorial_III.m - tutorial of the mcmc_simbio package
% demonstrating the estimation of parameters shared between two different
% circuits. In the language of the concurrent parameter inference
% problem (type 'help mcmc_info' into the command prompt window to 
% read more), we say that there are two
% network topologies, with each topology having one geometry associated
% with it.
%
% This example demonstrates a slightly more complex example of the
% concurrence feature of the mcmc_simbio Bayesian parameter inference
% toolbox. Here, we have two circuits: the constitutive gene expression
% circuit with model 
% 
% D + pol <-> D__pol  (kfdG, krdG)
% D__pol -> D + pol + protien v(kcp)
%
% and the tetR repression circuit with model
%
% D_T + P <-> D_T:P -> D_T + P + T (kfdT, krdT; kcp)
% D_G + P <-> D_G:P -> D_G + P + G (kfdG, krdG; kcp)
% 2 T <-> T2 (kfdimTet, krdimTet)
% D_G + T2 <-> D_G:T2 (kfseqTet, krseqTet)
%
% Here each model is a different topology with only one
% geometry associated with it. The paramters that are shared 
% between the two topologies are: kfdG, krdG, kcp, pol. The remaining
% parameters are specific to the topology-geometry pair they appear in (in
% this case all the remaining parameters appear in the tetR repression
% circuit topology). Furthermore, we set the forward rate parameters in
% all the reversible reaction to be fixed parameters, and therefore only
% estimate the reverse rate parameters. 

%% Initializing the toolbox
% If you have not already initialized the txtlsim and mcmc_simbio
% toolboxes, initialize them by running the txtl_init and mcmc_init
% commands in the command line. You need your working directory to be the
% main directory where the txtlsim toolbox is stored (i.e., the directory
% in which directories like core, components, mcmc_simbio etc are stored).

%% Run project_init 
% This creates a directory within the projects directory where 
% the results of the simulation will be stored. The name of the directory
% will be the same as the name of this file (proj_mcmc_tutorial_III, in this 
% case). it also creates a timestamped subdirectory within this directory, 
% where the actual results are stored. If the top level directory
% (proj_mcmc_tutorial_III) already exists, then only the subdirectory is
% created. 

[tstamp, projdir, st] = project_init;
prevtstamp = '20190201_170033'

 delete(gcp('nocreate'))
 parpool(46)
%% Define the MATLAB Simbiology model 
% We use the file model_protein3.m to define a constitutive gene expression 
% model using a single enzymatic step. 
m_constgfp = model_protein3;
m_tetRrep = model_tetR_repression1;

%% Defining the experiment / model arrangement. 
% We can define the experimental setup and how it is related to data, the
% Simbiology model and the estimation problem using what we call an
% mcmc_info struct. For this example, we will be using an
% mcmc_info_constgfp3tetR1.m file to generate the mcmc_info struct that we need
% to define our parameter inference problem. 
% Please enter 'help mcmc_info' into the command window prompt to read more
% about this struct. Also, open the mcmc_info_constgfp3tetR1 file (enter 'edit
% mcmc_info_constgfp3tetR1' into the command prompt) to learn how it is set up. 
mcmc_info = mcmc_info_constgfp3tetR1(m_constgfp, m_tetRrep);

%% Creating artificial data to fit the model to. 
% Instead of using real data, we will create artificial data for
% demonstration purposes. We will use the data_artificial_v2 fucntion to do
% this. 

% Get the model_info struct needed to generate the artificial data 
mi = mcmc_info.model_info; 

% A list of nominal parameter values to use to generate the data. 

cpol = 100; % nM
rkfdG = 5; % nM-1s-1
rkrdG = 300; % s-1
rkfdT = 5;
rkrdT = 300;
rkcp = 0.012; %s-1
rkfdimTet = 20; % nM-1s-1
rkrdimTet = 10; % s-1
rkfseqTet = 20; % nM-1s-1
rkrseqTet = 10; % s-1

% Arrange the parameters in a log transformed vector. 
masterVector = log([...
rkfdG 
rkrdG
rkfdT
rkrdT
rkfdimTet
rkrdimTet
rkfseqTet
rkrseqTet
rkcp
cpol]);

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

di = data_artificial_v2(...
    {m_constgfp, m_tetRrep},... % the two model objects
    {0:180:7200, 0:180:7200},... % time vectors for the two data sets
    {mi(1).measuredSpecies, mi(2).measuredSpecies},... 
    ...                 % measured species setup in mcmc_info.model_info
    {mi(1).dosedNames, mi(2).dosedNames},... % dosed species 
    {mi(1).dosedVals, mi(2).dosedVals,},... % dosing values
    {mi(1).namesUnord, mi(2).namesUnord},... 
    ...                 % names of species and parameters to set in each model
    {exp(masterVector([1 2 9 10])), exp(masterVector)}); 
                        % values to use for the names in namesUnord. 

 
da_constgfp = di(1).dataArray;
da_tetRrep = di(2).dataArray;
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

marray = mcmc_get_walkers({prevtstamp}, ...
    {ceil((SS.mcmc_info.runsim_info.nIter)/4):(SS.mcmc_info.runsim_info.nIter)},...
    projdir);
% assume the projdir where this data is stored is the same one as the
% one created at the start of this file

%%
pID = 1:length(mai.estNames);
marray_cut = mcmc_cut(marray, pID, flipud((mai.paramRanges)'));
if size(marray_cut, 2) < ri.nW
    warning('too few initial points, using a few timesteps from previous runs to create initial walker positions.');
    walker_timepoints = ceil(linspace(ceil(size(marray_cut,3))/4, size(marray_cut,3), ceil(ri.nW/size(marray_cut, 2))))
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

% mi = mcmc_runsim_v2(tstamp, projdir, di, mcmc_info,...
%    'InitialDistribution', 'LHS', 'multiplier', 2,...
%    'pausemode', true); 

% 'InitialDistribution', 'gaussian'
% 
%%  plot stuff 
% 
% These functions simply generate some standard plots from the data that is
% saved in the timestamped subdirectory of the directory specified in
% projdir. You can open that directory to view the results, including a log
% file. 
% 
% tstamptouse = tstamp; 
% marray = mcmc_get_walkers({tstamptouse}, {1:ri.nIter}, projdir);
% 
% % plot parameter distribution corner plot, and markov chains. 
% mcmc_plot(marray, mai.estNames,...
%     'savematlabfig', true, 'savejpeg', true,...
%     'projdir', projdir, 'tstamp', tstamptouse,...
%     'extrafignamestring', '_tutorialIII');
% 
% % plot individual trajectories of the data and the model fits for both
% % models. 
% titls = {'dG 10';'dG 30';'dG 60';};
% lgds = {};
% mvarray = masterVecArray(marray, mai);
% marrayOrd = mvarray(mi(1).paramMaps(mi(1).orderingIx, 1),:,:);
% fhandle = mcmc_trajectories(mi(1).emo, di(1), mi(1), marrayOrd,...
%     titls, lgds,...
%     'SimMode', 'curves', 'savematlabfig', true, 'savejpeg', true,...
%     'projdir', projdir, 'tstamp', tstamptouse, 'extrafignamestring',...
%     '_contgfp');
% marrayOrd = mvarray(mi(2).paramMaps(mi(2).orderingIx, 1),:,:);
% titls = {'dG 10 dT 0.1';'dG 30 dT 0.1';'dG 10 dT 2';'dG 30 dT 2';...
%     'dG 10 dT 8';'dG 30 dT 8';};
% fhandle = mcmc_trajectories(mi(2).emo, di(2), mi(2), marrayOrd,...
%     titls, lgds,...
%     'SimMode', 'curves', 'savematlabfig', true, 'savejpeg', true,...
%     'projdir', projdir, 'tstamp', tstamptouse, 'extrafignamestring',...
%     '_tetRrep');

%  Vipul Singhal, 
%  California Institute of Technology
%  2018
