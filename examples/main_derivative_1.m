clear all
close all
%% main_derivative_1
% first approach to implement derivative part:
% * plasR_ptet promoter activated by OC12HSL:protein lasR, repressed by tetR
% * OC12HSL:protein lasR is added continuously 
% * high DNA amount necessary
% * seems to work better with normal utr1 instead of (weak) utr14
% * Can we realize the ramp input in TXTL?
%% Negative Autoregulation - A simple example
% Here we demonstrate the setup of a genetic circuit where a transcription
% factor represses its own expression. 
% 
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
% ``E30VNPRL'' refers to a configuration file 
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');

u_slopes=[0,0.05,0.1,0.3,1,5];

for i = 1:length(u_slopes)
% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands, and all the relevant reactions
txtl_add_dna(tube3, 'plasR_ptet(50)', 'utr1(20)', 'tetR(1200))', 20, 'plasmid');
% txtl_add_dna(tube3, 'p70(50)', 'utr1(20)', 'deGFP(1000)', 1, 'plasmid');
% txtl_addspecies(tube3, 'OC12HSL:protein lasR', 10);
txtl_addreaction(tube3,' null -> [OC12HSL:protein lasR] ','MassAction',{'custom_prod_rate',u_slopes(i),'F'});
% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!

cs = getconfigset(Mobj);
set(cs.RuntimeOptions, 'StatesToLog', 'all');
tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;

%% plot the result
% The following function plots the proteins, RNA and resources in the
% toolbox. In the next section we delve deeper into the object oriented
% structure of the model, and how to plot arbitrary species in the model. 
% txtl_plot(simData,Mobj);

%% Model Structure
% The model is organized as a model object, with sub objects specifying
% Parameters, Reactions, Species, etc. Type in 
% Mobj

%%
% There is one comaprtment, 2 events, 73 parameters, 47 Reactions, no rules
% and 43 Species in the toolbox. We can explore further by typing,
% for example, 
% Mobj.Species

%%
% We see that there are 43 species in the model, and they have somewhat 
% different syntax for specification. Proteins, RNA and DNA generally follow the convention
% protein CDS, RNA 5'UTR--CDS, DNA promoter--5' UTR--CDS, with
% variations possible. There are also simply named `core' species like
% RNAP, Ribo, RNase, etc. Finally we denote bound complexes with a colon,
% for example, Species 1:Species 2. 
%% 
% We also see that each of them has certain
% other associated properties. You can explore further by accessing
% individual species using their index, and using the `get' and `set' commands
% to get and set the properties of the species. For example, try typing 
% Mobj.Species(1)

%% 
% This gives you the first species in the model. You can find out what
% properties as associated with this species by typing in 
% get(Mobj.Species(1))

%%
% and then using the set command to set its initial concentration to 50
% units:
% set(Mobj.Species(1), 'InitialAmount', 50)

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

Inducerindex = findspecies(Mobj, 'OC12HSL:protein lasR','withInComplex');
reporterindex = findspecies(Mobj, 'RNA utr1--tetR','withInComplex');
intermediateindex = findspecies(Mobj, 'protein tetR','withInComplex'); 
timeset{i}=simData.Time/3600;
uset{i}=sum(simData.data(:,Inducerindex),2);
yset{i}= sum(simData.data(:,reporterindex),2);
zset{i}= sum(simData.data(:,intermediateindex),2);


end

figure
subplot(1,3,1)
hold all
title('OC12HSL:protein lasR as input u')
ylabel('concentration, AU')
xlabel('time, AU')
for i = 1:length(uset)
    plot(timeset{i},uset{i});
end


subplot(1,3,2)
hold all
title('RNA utr1--tetR as output y')
ylabel('concentration, AU')
xlabel('time, AU')
for i = 1:length(uset)
    plot(timeset{i},yset{i});
end
legend(cellstr(num2str(u_slopes', '%.2f')))

subplot(1,3,3)
hold all
title(' protein tetR as output z')
ylabel('concentration, AU')
xlabel('time, AU')
for i = 1:length(uset)
    plot(timeset{i},zset{i});
end

% % 



%%

%% EMACS editor support (ignore)
%%
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
