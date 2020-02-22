%
function [simData,Mobj] = phosphor_insulator_function(extract,dna_sigma54,dna_NRII_L16R,dna_NRII_H139N,dna_NRI,dna_deGFP,tspan)

tube1 = txtl_extract(extract);
tube2 = txtl_buffer(extract);

tube3 = txtl_newtube('gene_expression');

% sigma factor
txtl_add_dna(tube3, ...
  'p70(50)', 'utr1(20)', 'sigma54(1434)', ...	% promoter, utr1, gene
   dna_sigma54, ...					% concentration (nM)
  'linear');					% type

% NRII_L16R KINASE
txtl_add_dna(tube3, ...
  'p70(50)', 'utr1(20)', 'NRII_L16R(1050)', ...	% promoter, utr1, gene
   dna_NRII_L16R, ...					% concentration (nM)
  'linear');					% type

% NRII_H139N
txtl_add_dna(tube3, ...
  'p70(50)', 'utr1(20)', 'NRII_H139N(1050)', ...	% promoter, utr1, gene
   dna_NRII_H139N, ...					% concentration (nM)
  'linear');					% type

% NRI
txtl_add_dna(tube3, ...
  'p70(50)', 'utr1(20)', 'NRI(1410)', ...	% promoter, utr1, gene
  dna_NRI, ...					% concentration (nM)
  'linear');					% type

txtl_add_dna(tube3, ...
  'pGlnA(308)', 'utr1(20)', 'deGFP(1000)', ...	% promoter, utr1, gene
   dna_deGFP, ...					% concentration (nM)
  'linear');					% type



Mobj = txtl_combine([tube1, tube2, tube3]);


tic
[simData] = txtl_runsim(Mobj,tspan*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;




