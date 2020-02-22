% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
% 
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
% 
close all 
clear all

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('ptet_charac');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'tetR(1000)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type
dna_deGFP = txtl_add_dna(tube3, ...
  'ptet(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
   2, ...					% concentration (nM)
  'plasmid');

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);
   txtl_addspecies(Mobj, 'aTc', 100);

cs = getconfigset(Mobj);
set(cs.RuntimeOptions, 'StatesToLog', 'all');
tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;



%% plot the result

 txtl_plot(simData,Mobj);
 

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
