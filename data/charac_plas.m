                               
close all 
% clear all
tube1 = txtl_extract('E30');
tube2 = txtl_buffer('E30');
tube3 = txtl_newtube('plas_charac');
dna_tetR = txtl_add_dna(tube3,  'placI(50)', 'rbs(20)', 'lasR(1000)',  1, 'plasmid'); % 163 in zach's index	
dna_deGFP = txtl_add_dna(tube3,  'plasR(50)', 'SNP10UTR1(20)', 'deGFP(1000)', 2,'plasmid'); % 196
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'OC12HSL', .1);
cs = getconfigset(Mobj);
set(cs.RuntimeOptions, 'StatesToLog', 'all');
tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;

%% plot the result

 txtl_plot(simData,Mobj);
 

 
cellOfSpecies1 = { 'AGTP','CUTP','AGTP_UNUSE'
    'DNA placI--rbs--lasR', 'RNAP70:DNA placI--rbs--lasR','term_RNAP70:DNA placI--rbs--lasR'
    'AGTP:RNAP70:DNA placI--rbs--lasR','CUTP:RNAP70:DNA placI--rbs--lasR','CUTP:AGTP:RNAP70:DNA placI--rbs--lasR'
    'RNAP70', 'RNAP28', 'RNAP'  };
 plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies1)  

            
cellOfSpecies2 = {'Ribo','AA','RNase'
                 'RNA rbs--lasR','Ribo:RNA rbs--lasR','AA:AGTP:Ribo:RNA rbs--lasR'
                 'RNA rbs--lasR:RNase','Ribo:RNA rbs--lasR:RNase', 'AA:AGTP:Ribo:RNA rbs--lasR:RNase'};
plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies2)


cellOfSpecies1 = {'DNA plasR--SNP10UTR1--deGFP', 'RNAP70:DNA plasR--SNP10UTR1--deGFP','term_RNAP70:DNA plasR--SNP10UTR1--deGFP:OC12HSL:protein lasR'
    'AGTP:RNAP70:DNA plasR--SNP10UTR1--deGFP:OC12HSL:protein lasR','CUTP:RNAP70:DNA plasR--SNP10UTR1--deGFP:OC12HSL:protein lasR','CUTP:AGTP:RNAP70:DNA plasR--SNP10UTR1--deGFP:OC12HSL:protein lasR'};
 plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies1)  

            
cellOfSpecies2 = {'OC12HSL', 'OC12HSL:protein lasR', 'RNAP70:DNA plasR--SNP10UTR1--deGFP:OC12HSL:protein lasR'
                  'RNA SNP10UTR1--deGFP','Ribo:RNA SNP10UTR1--deGFP','AA:AGTP:Ribo:RNA SNP10UTR1--deGFP'
                 'RNA SNP10UTR1--deGFP:RNase','Ribo:RNA SNP10UTR1--deGFP:RNase', 'AA:AGTP:Ribo:RNA SNP10UTR1--deGFP:RNase'};
plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies2)

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
