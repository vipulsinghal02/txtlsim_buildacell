%% Tutorial II
% tutorial_ii.m - Slightly more involved set of examples involving txtlsim. We
% recommend you read through the
% <https://vipulsinghal02.github.io/txtlsim_buildacell/tutorial.html
% first> tutorial. 
%
% Vipul Singhal, 18 Mar 2018
%
% In this tutorial, we explore and demonstrate some of the finer points of
% the txtlsim toolbox. 


%% Initializing the toolbox
% Remember to set the working directory to the trunk directory of the
% toolbox. 
% 
% Use this command to add the subdirectories needed to your matlab path. To
% be run each time you begin a new TXTL toolbox session. 
txtl_init;

%% IFFL example
% The code below can be used to set up an IFFL example with lasR as the
% activator, tetR as the repressor, and deGFP as the reporter. 
% 
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
% ``E2'' refers to a configuration file 
tube1 = txtl_extract('E3');
tube2 = txtl_buffer('E3');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('lastetIFFL');


% Define the DNA strands, and all the relevant reactions
txtl_add_dna(tube3, ...
  'plac(50)', 'utr1b(20)', 'lasR(1000)', 1, 'plasmid');	
txtl_add_dna(tube3, ...
  'plas(50)', 'utr1b(20)', 'tetR(1000)', 0.1, 'plasmid');	
txtl_add_dna(tube3, ...
  'plas_ptet(50)', 'utr1b(20)', 'deGFP(1000)', 1, 'plasmid');	

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

txtl_addspecies(Mobj, 'OC12HSL', 1000);
txtl_addspecies(Mobj, 'aTc', 100);

% Run a simulaton
tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;

%% Plot the result
% The following function plots the proteins, RNA and resources in the
% toolbox. In the next section we delve deeper into the object oriented
% structure of the model, and how to plot arbitrary species in the model. 
txtl_plot(simData,Mobj);

%% Model Structure
% The model is organized as a model object, with sub objects specifying
% Parameters, Reactions, Species, etc. Type in 
Mobj

%% 
% and the individual properties of the model object may be explored by
% typing, for example, 
Mobj.Events


%% The effect of varying aTc
% we will vary the aTc inducer in the IFFL circuit, and plot the
% results of these variations. We will use MATLAB cell arrays to store the
% variations of the models. 
close all
initATC = [1 10 100 1000];
m = cell(length(initATC), 1);
simData = cell(length(initATC), 1);
t_ode = cell(length(initATC), 1);
x_ode = cell(length(initATC), 1);

tube1 = txtl_extract('E3');
tube2 = txtl_buffer('E3');
for i = 1:length(initATC)
    tube3 = txtl_newtube('lastetIFFL');
    txtl_add_dna(tube3, ...
      'plac(50)', 'utr1b(20)', 'lasR(1000)', 1, 'plasmid');	
    txtl_add_dna(tube3, ...
      'plas(50)', 'utr1b(20)', 'tetR(1000)', 0.1, 'plasmid');	
    txtl_add_dna(tube3, ...
      'plas_ptet(50)', 'utr1b(20)', 'deGFP(1000)', 1, 'plasmid');	
    m{i} = txtl_combine([tube1, tube2, tube3]);
    txtl_addspecies(m{i}, 'OC12HSL', 1000);
    txtl_addspecies(m{i}, 'aTc', initATC(i));
    [simData{i}] = txtl_runsim(m{i},6*60*60);
    t_ode{i} = simData{i}.Time;
    x_ode{i} = simData{i}.Data;
end


speciesList = {'aTc'
'protein tetRdimer' 
'RNA utr1b--deGFP'
'protein deGFP*'};

plotCustomSpecies2(m, x_ode,t_ode, speciesList, {'1nM', '10nM', '100nM', '1000nM'})
