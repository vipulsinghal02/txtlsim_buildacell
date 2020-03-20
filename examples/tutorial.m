%% Your very first tutorial on the txtlsim toolbox
% tutorial.m - basic usage of the TXTL modeling toolbox
%
% Vipul Singhal, 28 July 2017
%
% This file contains a simple tutorial of the TXTL modeling toolbox. You
% will learn about setting up a constitutive expression circuit, and a
% negative autoregulation circuit. You will also learn the basics of
% plotting the results, creating simple variations of a circuit, and of how
% the model is structured. 

%% Initializing the toolbox
% Remember to set the working directory to the trunk directory of the
% toolbox. The trunk directory is where folders like ``core'', and
% ``components''live. Here is a snapshot of this directory on the author's
% computer. 
% 
% <<trunk_dir.png>>
% 
% !TODO test this. does the image show up online? Probably not. 
% 
% Use this command to add the subdirectories needed to your matlab path. To
% be run each time you begin a new TXTL toolbox session. 
txtl_init;

%% Negative Autoregulation - A simple example
% Here we demonstrate the setup of a genetic circuit where a transcription
% factor represses its own expression. 
% 
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
% ``E2'' refers to a configuration file 
tube1 = txtl_extract('E2');
tube2 = txtl_buffer('E2');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands, and all the relevant reactions
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)', 'tetR(1200)', 1, 'plasmid');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)', 'deGFP(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% Run a simulaton
%   
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!

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

%% Model structure
% The model is organized as a model object, with sub objects specifying
% Parameters, Reactions, Species, etc. Type in 
Mobj

%%
% We can see the number of instances of the various subclasses of the 
% model object. We can explore further by typing 
Mobj.Species

%%
% Proteins, RNA and DNA generally follow the naming convention
% descirbed in the paper, with, for example, DNA specified with ``DNA
% promoter--5' UTR--CDS''.
%
% There are also simply named `core' species like
% RNAP, Ribo, RNase, etc. Finally we denote bind together speices into 
% complexes with a colon, for example, Species 1:Species 2.  
%% 
% We also see that each of them has certain
% other associated properties. You can explore further by accessing
% individual species using their index, and using the `get' and `set' commands
% to get and set the properties of the species. For example, try typing 
Mobj.Species(1)

%% 
% This gives you the first species in the model. You can find out what
% properties as associated with this species by typing in 
get(Mobj.Species(1))

%%
% and then using the set command to set its initial concentration to 50
% units:
set(Mobj.Species(1), 'InitialAmount', 50)

%%
% Learn more about the get and set commands by typing in 
%%
%   help get
%   help set

%% 
% You may read more about how model objects are arranged in Simbiology by
% working through the
% <https://www.mathworks.com/help/simbio/gs/simbiology-command-line-tutorial.html
% Tutorial>. Feel free to browse the reactions and other subproperties by
% individually typing in commands shown below: 
%% 
% The reactions of the model object: 
  Mobj.reactions
%% 
% The properties of the reactions that can be queried:
  get(Mobj.Reactions)
%% 
% The properties of the first reaction:
  get(Mobj.Reactions(1))
%% 
% The reaction rate of the first reaction:
  Mobj.Reactions(1).ReactionRate
%% 
% The Kinetic law object used by the first reaction:
  Mobj.Reactions(1).KineticLaw
%% 
% All the visible properties of the kinetic law object of the first
% reaction:
  get(Mobj.Reactions(1).KineticLaw)
%%
% and so on.
%% Plotting individual species
% You can also plot the trajectories of any of the species in the model.
% Use the function findspecies to get the index of the species object of interest. For
% example, if you want to plot the trajectory of the dimerized tetR
% protein, you could type in
tetRindex = findspecies(Mobj, 'protein tetRdimer');
figure
plot(simData.Time/3600, simData.data(:,tetRindex));
title('Dimerized tetR concentration')
ylabel('concentration, nM')
xlabel('time, h')
curraxis = axis; 
axis([curraxis(1:2) 0 curraxis(4)])

