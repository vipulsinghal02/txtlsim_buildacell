function [simData,Mobj] = incoherent_ff_loop_function(extract,p70_AraC,pBAD_tetR,pBAD_ptet_deGFP,tspan,protein_ClpX,arabinose,aTc)


tube1 = txtl_extract(extract);
tube2 = txtl_buffer(extract);

tube3 = txtl_newtube('circuit_closed_loop_withClpX');

txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC(600)',p70_AraC, 'plasmid');
txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(600)', pBAD_tetR, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)-lva(20)',pBAD_ptet_deGFP, 'plasmid');

txtl_addspecies(tube3, 'arabinose', arabinose);
txtl_addspecies(tube3, 'aTc', aTc);
txtl_addspecies(tube3, 'protein ClpX', protein_ClpX);
% set up well_a1
Mobj = txtl_combine([tube1, tube2, tube3]);

%% Run a simulation
 simulationTime = tspan*60*60;
tic
simData = txtl_runsim(Mobj,simulationTime);
toc




end