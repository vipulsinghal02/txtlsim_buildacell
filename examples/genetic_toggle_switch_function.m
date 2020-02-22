
%
% This file contains a simple example of setting up a TXTL simulation
% for a genetic switch. It shows the bistability of the circuit, and the
% switching that can be accomplished by the tetR repressing inducer aTc and
% the lacI repressing inducer IPTG. 

% the dimerization and tetramerization rates of lacI and tetR here are
% high. need to check what the actual values are, since that is intimitaly
% related to the working of this switch. 
% another thing to be played with (and hence should be told to anyone who
% uses this simulation, is that TX and TL rates might be mo

function [simData,Mobj] = genetic_toggle_switch_function(extract,ptet_DNA,placI_DNA,tspan,IPTG_amount,atc_amount)


% set up the tubes
tube1 = txtl_extract(extract);
tube2 = txtl_buffer(extract);
tube3 = txtl_newtube('genetics_switch');

% add dna to tube3
dna_lacI = txtl_add_dna(tube3,'ptet(50)', 'rbs(20)', 'lacI(647)', ptet_DNA, 'plasmid');
dna_tetR = txtl_add_dna(tube3, 'plac(50)', 'rbs(20)', 'tetR(647)', placI_DNA, 'plasmid');
% dna_deGFP = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'deGFP(1000)', ptet_DNA, 'plasmid');
% dna_deCFP = txtl_add_dna(tube3, 'placI(50)', 'rbs(20)', 'deCFP(1000)', placI_DNA, 'plasmid');

% combine extact, buffer and dna into one tube
Mobj = txtl_combine([tube1, tube2, tube3]);

% add inducers
txtl_addspecies(Mobj, 'IPTG',IPTG_amount);
txtl_addspecies(Mobj, 'aTc',atc_amount);

% set up simulation
simulationTime = tspan*60*60; % hours

% simulate the reaction
tic 
simData = txtl_runsim(Mobj,simulationTime);
toc

end
