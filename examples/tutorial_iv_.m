%% Tutorial IV
% tutorial_iv.m - Circuit components and variations. 
% 
% In this tutorial, we go over some of the various circuit components that
% can be used to build circuit models. We also discuss variations like
% using linear DNA and protection from DNA degradation, and ClpX. 
% We also discuss the parameters and reactions in the
% models, so that you may gain greater familiarity with the 
% 
% 
% We  recommend you read through tutorials i and ii before trying 
% this one out.
% All the tutorials can be found on the index page,
% <https://vipulsinghal02.github.io/txtlsim_buildacell/
% here>.
% 
% Vipul Singhal

%% todos
% # linear DNA degradation example. -- recBCD, gamS. Zach paper. 
% # ClpX degradation example. (from tech report)
% # activator example + inducer (use our characterization)
% # combinatorial promoter example + inducers. (use our characterization)
% #repressor example (use our characterization, and lacI too, find a characterization of lacI in TXTL/ or just relative to tetR.)
% # rbs varying example (use the mutalik paper to estimate the RBSs)
% # constitutive promoter example (use the anderson promoter library to
% # estimate promoter strength )
% # genetic toggle with the data collected. (anushka, enrique)
% # Can look at the TXTL2.0 paper for possible data.
% maybe wan to split this into multiple smaller links? No need. just have a
% contents page. 
% 
 


%% Initializing the toolbox
% Remember to set the working directory to the trunk directory of the
% toolbox. 
% 
% Use this command to add the subdirectories needed to your matlab path. To
% be run each time you begin a new TXTL toolbox session. 
txtl_init;

%% IFFL example
% The code below can be used to set up the constitutive GFP production
% example. 
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
  'plac(50)', 'utr1(20)', 'lasR(1000)', 5, 'plasmid');	
txtl_add_dna(tube3, ...
  'plas(50)', 'utr1(20)', 'tetR(1000)', 0.2, 'plasmid');	
txtl_add_dna(tube3, ...
  'plas_ptet(50)', 'utr1(20)', 'deGFP(1000)', 2, 'plasmid');	
m = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(m, 'OC12HSL', 1000);
txtl_addspecies(m, 'aTc', 1000);

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% Run a simulaton
% 

tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;

%% plot the result
% The following function plots the proteins, RNA and resources in the
% toolbox. In the next section we delve deeper into the object oriented
% structure of the model, and how to plot arbitrary species in the model. 
txtl_plot(simData,Mobj);

%% Model Structure
% The model is organized as a model object, with sub objects specifying
% Parameters, Reactions, Species, etc. Type in 
Mobj

%%
% We can see the number of instances of the various subclasses of the 
% model object. We can explore further by typing 
Mobj.Species

%%
% Proteins, RNA and DNA generally follow the convention
% protein CDS, RNA 5'UTR--CDS, DNA promoter--5' UTR--CDS, with
% variations possible. There are also simply named `core' species like
% RNAP, Ribo, RNase, etc. Finally we denote bound complexes with a colon,
% for example, Species 1:Species 2. 
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
% individually typing in commands like 
%%
% 
%   Mobj.reactions
%   get(Mobj.Reactions)
%   get(Mobj.Reactions(1))
%   Mobj.Reactions(1).ReactionRate
%   Mobj.Reactions(1).KineticLaw
%   get(Mobj.Reactions(1).KineticLaw)
%%
% and so on.
%% Plotting individual species
% You can also plot the trajectories of any of the species in the model.
% Use the function findspecies to get the index of the species object of interest. For
% example, if you want to plot the trajectory of the dimerized tetR
% protein, you could type in

tetRindex = findspecies(Mobj, 'protein deGFP');
figure
plot(simData.Time/3600, simData.data(:,tetRindex));
title('Un-matured protein concentration')
ylabel('concentration, nM')
xlabel('time, h')
curraxis = axis; 
axis([curraxis(1:2) 0 curraxis(4)])

