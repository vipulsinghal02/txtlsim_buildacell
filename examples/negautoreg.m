% negautoreg.m - negative autoregulation example
% R. M. Murray, 8 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a negatively autoregulated gene.  The constants for this example
% come from the Simbiology toolbox example page:
%
%    http://www.mathworks.com/help/toolbox/simbio/gs/fp58748.html
%
tube1 = txtl_extract('E2');
tube2 = txtl_buffer('E2');
tube3 = txtl_newtube('negautoreg');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'tetR(1200)', 10, 'plasmid');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'deGFP(1000)', 10, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'aTc', 1000);
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
