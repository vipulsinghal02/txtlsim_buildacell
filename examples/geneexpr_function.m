
function [Mobj,simData] =  geneexpr_function(extract,promoter,dna_amount,tspan,varargin)
% geneexpr_function  simulates gene expression under different conditions.
%   [Mobj,simData] = geneexpr_function(extract,promoter,dna_amount,tspan)
%   extract states the extract config file from the config file directory,
%   promoter specify the promoter and dna_amount sets DNA concetration in
%   the simulation. tspan tells to the ode solver, the simulation time (in hours!)
%
%   See also geneexpr


% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract(extract);
tube2 = txtl_buffer(extract);

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_add_dna(tube3, ...
  promoter, 'utr1(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
 dna_amount, ...					% concentration (nM)
  'plasmid');					% type

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% Run a simulation
simulationTime = tspan*60*60;


%run
if nargin > 4
    
    % call runsim to set up reactions with simulation
    txtl_runsim(Mobj);
     
    % read out the desired parameter values and set them up
    %reactionindex, paramter num, value
    paramLoc = varargin{1}; 
    % that could be done with cellfun, but in this way it is easier to
    % undertand 
    for k = 1:size(paramLoc,1)
        set(Mobj.Reactions(paramLoc{k,1}).KineticLaw.Parameters(paramLoc{k,2}),'Value',paramLoc{k,3})
    end
%     sbioaccelerate(Mobj, configsetObj)
    [simData] = txtl_runsim(Mobj,simulationTime);
    
else
    % normal operation mode
    [simData] = txtl_runsim(Mobj,simulationTime);
end



end
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:

