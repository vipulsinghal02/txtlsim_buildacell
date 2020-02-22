
%%
% *negative auto regulation function form*
% input:
% * extract: name of the extract
% * dna_amount (plasmid): in the reaction
% * tspan: the simulation time in *hours*
% * atc_amount: amount of ATC induder in the reaction
% 
% output:
% * simulation timepoint
% * simulation data
% * Simbiology model obj

function [t_ode, x_ode, Mobj] =  negautoreg_function(extract,dna_amount,tspan,atc_amount)

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract(extract);
tube2 = txtl_buffer(extract);

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('negautoreg');


% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)', 'tetR(1200)', dna_amount, 'plasmid');



%
% Next we have to set up the reactions that describe how the circuit
% works.  Transcription and translation are already included above, so
% we just need to include protein-protein and protein-DNA interactions.
%
% Note that the commands in this section are standard Simbiology commands,
% so you can put anything you want here.
%

% No additional reactions required for this circuit
% tetR-DNA interactions are automatically included in tetR setup

%
% Describe the actual experiment that we want to run.  This includes 
% combining the various tubes and also adding any additional inducers
% or purified proteins that you want to include in the run.
%

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);
if ~isempty(atc_amount)
 txtl_addspecies(Mobj, 'aTc', atc_amount);
end
 

% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
simulationTime = tspan*60*60;
[t_ode, x_ode] = txtl_runsim(Mobj, simulationTime);


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
