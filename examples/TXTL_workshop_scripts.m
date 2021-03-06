% TXTL modeling workshop examples
% Vipul Singhal
% California Institute of Technology

%% Geneexpr without protein degradation
clear all; close all; clc
tube1 = txtl_extract('E2');
tube2 = txtl_buffer('E2');
tube3 = txtl_newtube('geneexpr_noClpX');
txtl_add_dna(tube3, 'p70(50)', 'utr1(20)',...
    'deGFP(1000)-lva(20)', 10, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);

%% Geneexpr with protein degradation
clear all;
tube1 = txtl_extract('E2');
tube2 = txtl_buffer('E2');
tube3 = txtl_newtube('geneexpr_ClpX');
txtl_add_dna(tube3, 'p70(50)', 'utr1(20)',...
    'deGFP-lva(1000)', 10, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'protein ClpX*', 100);
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);

%% Geneexpr without protein degradation, but with linea DNA (degraded by RecBCD) 
tube1 = txtl_extract('E2');
tube2 = txtl_buffer('E2');
tube3 = txtl_newtube('geneexpr_noClpX');
txtl_add_dna(tube3, 'p70(50)', 'utr1(20)',...
    'deGFP(1000)-lva(20)', 10, 'linear');
Mobj = txtl_combine([tube1, tube2, tube3]);
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);

%% Geneexpr without protein degradation, but with linea DNA (degraded by RecBCD) 
% but protected by gamS

tube1 = txtl_extract('E2');
tube2 = txtl_buffer('E2');
tube3 = txtl_newtube('geneexpr_noClpX');
txtl_add_dna(tube3, 'p70(50)', 'utr1(20)',...
    'deGFP(1000)-lva(20)', 10, 'linear');
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'protein gamS', 3500);
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);
%% Negative Autoreg with plasmid DNA + aTc induceer
clear all; close all; clc
% Plasmid DNA
tube1 = txtl_extract('E2');
tube2 = txtl_buffer('E2');
tube3 = txtl_newtube('negautoreg_plasmid');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'tetR(1200)', 10, 'plasmid');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'deGFP(1000)', 10, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'aTc', 50); % try changing the aTc value to 0, 5, 50, 500
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);

%% linear unprotected DNA
tube3 = txtl_newtube('negautoreg_linear_noGamS');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'tetR(1200)', 10, 'linear');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'deGFP(1000)', 10, 'linear');
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'aTc', 50);
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);

%% linear protected DNA
tube3 = txtl_newtube('negautoreg_linear_GamS');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'tetR(1200)', 10, 'linear');
txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
    'deGFP(1000)', 10, 'linear');
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'aTc', 50);
txtl_addspecies(Mobj, 'protein gamS', 12000);
[simData] = txtl_runsim(Mobj,14*60*60);
txtl_plot(simData,Mobj);
%% !TODO:
% %% Inducer Runs: linear DNA
% 
% % initialize arrays
% tetRconc = [1 1 1 5 5 5];
% aTcconc = [10 500 3000 10 500 3000];
% t_cellarray = cell(6,1);
% x_cellarray = cell(6,1);    
% Mobj_cellarray = cell(6,1);
% 
% % run code over all initial conditions, saving data in cell arrays
% for i = 1:6
% tube1 = txtl_extract('E30VNPRL');
% tube2 = txtl_buffer('E30VNPRL');
% tube3 = txtl_newtube('induction');
% txtl_add_dna(tube3, 'plac(50)', 'utr1(20)',...
%     'tetR(1200)', tetRconc(i), 'linear');
% txtl_add_dna(tube3, 'ptet(50)', 'utr1(20)',...
%     'deGFP(1000)', 10, 'linear');
% Mobj_cellarray{i} = txtl_combine([tube1, tube2, tube3]);
% txtl_addspecies(Mobj_cellarray{i}, 'aTc', aTcconc(i));
% txtl_addspecies(Mobj_cellarray{i}, 'protein gamS', 3500);
% [t_cellarray{i}, x_cellarray{i}] = txtl_runsim(Mobj_cellarray{i},14*60*60);
% end
% 
% % plot species of interest
% speciesToPlot = {'protein deGFP*', 'protein deGFP'
%     'protein tetRdimer', '2 aTc:protein tetRdimer', };
% legendList = {'tetR DNA = 1nM, aTc = 10nM',...
%     'tetR DNA = 1nM, aTc = 500nM',...
%     'tetR DNA = 1nM, aTc = 3000nM',...
%     'tetR DNA = 5nM, aTc = 10nM',...
%     'tetR DNA = 5nM, aTc = 5000nM',...
%     'tetR DNA = 5nM, aTc = 3000nM'};
% plotCustomSpecies2(Mobj_cellarray, x_cellarray, t_cellarray, speciesToPlot, legendList);
% 
% %% IFFL runs: plamid DNA
% 
% 
% tube1 = txtl_extract('E30VNPRL');
% tube2 = txtl_buffer('E30VNPRL');
% tube3 = txtl_newtube('IFFL_Sguo_plasmid');
% 
% dna_lasR = txtl_add_dna(tube3,  'plac(50)', 'utr1(20)', 'lasR(1000)',  1, 'plasmid'); % 163 in zach's index	
% dna_tetR = txtl_add_dna(tube3,  'plas(50)', 'utr1(20)', 'tetR(1000)', 0.1,'plasmid'); % 196
% dna_deGFP = txtl_add_dna(tube3,  'plas_ptet(50)', 'utr1(20)', 'deGFP(1000)-lva', 2,'plasmid');
% 
% Mobj = txtl_combine([tube1, tube2, tube3]);
% % txtl_addspecies(Mobj, 'protein gamS', 3500);
% txtl_addspecies(Mobj, 'OC12HSL', 2000);
% txtl_addspecies(Mobj, 'aTc', 500); %try 10, 500 and 5000
% % txtl_addspecies(Mobj, 'arabinose', 1000);
% txtl_addspecies(Mobj, 'protein ClpX*', 100);
% tic
% [simData] = txtl_runsim(Mobj,10*60*60);
% 
% toc
% txtl_plot(simData,Mobj);
% speciesToPlot = {'protein deGFP-lva*', 'protein deGFP-lva', 'protein deGFP-lva***'
%     'protein tetRdimer', '2 aTc:protein tetRdimer','aTc'
%     'protein lasR', 'OC12HSL:protein lasR', 'OC12HSL'};
% t_ode = simData.Time;
% x_ode = simData.Data;
% plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, speciesToPlot);
% 
% 
% %% IFFL with linear DNA
% tube1 = txtl_extract('E30VNPRL');
% tube2 = txtl_buffer('E30VNPRL');
% tube3 = txtl_newtube('IFFL_Sguo_linear');
% 
% dna_lasR = txtl_add_dna(tube3,  'plac(50)', 'utr1(20)', 'lasR(1000)',  10, 'linear'); % 163 in zach's index	
% dna_tetR = txtl_add_dna(tube3,  'plas(50)', 'utr1(20)', 'tetR(1200)', 10,'linear'); % 196
% dna_deGFP = txtl_add_dna(tube3,  'plas_ptet(50)', 'utr1(20)', 'deGFP(1000)-lva', 10,'linear');
% 
% Mobj = txtl_combine([tube1, tube2, tube3]);
% txtl_addspecies(Mobj, 'protein gamS', 100000);
% txtl_addspecies(Mobj, 'OC12HSL', 1000);
% txtl_addspecies(Mobj, 'aTc', 1000);
% txtl_addspecies(Mobj, 'protein ClpX*', 100);
% tic
% [simData] = txtl_runsim(Mobj,10*60*60);
% toc
% txtl_plot(simData,Mobj);
% speciesToPlot = {'protein deGFP-lva*', 'protein deGFP-lva', 'protein deGFP-lva***'
%     'protein tetRdimer', '2 aTc:protein tetRdimer','aTc'
%     'protein lasR', 'OC12HSL:protein lasR', 'OC12HSL'};
% t_ode = simData.Time;
% x_ode = simData.Data;
% plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, speciesToPlot);