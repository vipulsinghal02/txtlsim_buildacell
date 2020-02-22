
% IFFL as tested in lab by S. Guo
% VS 2013

%% clean up


close all

%% no clpX
% Set up the standard TXTL tubes


tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');

tube3 = txtl_newtube('circuit_closed_loop_withClpX');

txtl_add_dna(tube3, 'p70(50)', 'utr1(20)', 'AraC(600)',0.5*4.5, 'plasmid');
txtl_add_dna(tube3, 'pBAD(50)', 'utr1(20)', 'tetR(600)', 2*4.5, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'utr1(20)', 'deGFP(1000)-lva(20)',4*4.5, 'plasmid');
txtl_add_dna(tube3,'p70(50)', 'utr1(20)', 'ClpX(600)',0.1*4.5, 'plasmid');

txtl_addspecies(tube3, 'arabinose', 10000);
txtl_addspecies(tube3, 'aTc', 1000);
ClpXToAdd = 80;
% set up well_a1
Mobj = txtl_combine([tube1, tube2, tube3]);

%% Run a simulation
 simulationTime = 0.01*60*60;

simData = txtl_runsim(Mobj,simulationTime);

t_ode = simData.Time;
x_ode = simData.Data;

iClpX = findspecies(Mobj, 'protein ClpX');
x_ode(end,iClpX) = x_ode(end,iClpX)+ClpXToAdd;%

 simulationTime = 8*60*60;
 %!TODO: vs 06/18/14 found a bug in runsim here. line 224. look into this. 
 [t_ode,x_ode] = txtl_runsim(Mobj,simulationTime, t_ode,x_ode);

defaultGroups = txtl_getDefaultPlotDataStruct();
defaultGroups(2).SpeciesToPlot = {'protein deGFP-lva*'};

txtl_plot(t_ode,x_ode,Mobj,defaultGroups);



%%
% % % % 
% % % % 
% % % % cellOfSpecies = {'RNAP70:DNA pBAD_ptet--utr1--deGFP-lva:protein tetRdimer', 'RNAP70:DNA pBAD_ptet--utr1--deGFP-lva:protein tetRdimer:arabinose:protein AraC'
% % % %                  'NTP:RNAP70:DNA pBAD_ptet--utr1--deGFP-lva:protein tetRdimer','NTP:RNAP70:DNA pBAD_ptet--utr1--deGFP-lva:protein tetRdimer:arabinose:protein AraC'
% % % %                  'RNAP70:DNA pBAD_ptet--utr1--deGFP-lva:arabinose:protein AraC', 'RNAP70:DNA pBAD--utr1--tetR:arabinose:protein AraC'
% % % %                  'NTP:RNAP70:DNA pBAD_ptet--utr1--deGFP-lva:arabinose:protein AraC', 'NTP:RNAP70:DNA pBAD--utr1--tetR:arabinose:protein AraC'};
% % % % plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies)
% % % % 
% % % % cellOfSpecies = {'protein deGFP-lva*', 'protein tetR','protein ClpX*'
% % % %                  'protein AraC','protein tetRdimer','protein deGFP-lva*:protein ClpX*'
% % % %                  'arabinose', 'aTc','protein deGFP-lva***'
% % % %                  'arabinose:protein AraC', '2 aTc:protein tetRdimer','protein ClpX'};
% % % % plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies)
% % % % 
% % % % cellOfSpecies = {'RNAP', 'protein sigma70','Ribo'
% % % %                  'RNAP70','RNase','NTP'
% % % %                  'ATP', 'AA','protein deGFP-lva'};
% % % % plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies)


%% Redo, but with varying levels of aTc

%initialize arrays
clear all

Param = 'tetR';
vParam = [0 0.1 0.5 2 5]; 
nParam = length(vParam);
Mobj = cell(nParam,1);
simData = cell(nParam,1);
t_ode = cell(nParam,1);
x_ode = cell(nParam,1);
iGFP = cell(nParam,1);
iAraC = cell(nParam,1);
itetR = cell(nParam,1);
iATC = cell(nParam,1);


%Run code
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');



for i = 1:nParam
tube3 = txtl_newtube('circuit_closed_loop_withClpX');
txtl_add_dna(tube3, 'p70(50)', 'utr1(20)', 'AraC(600)',0.5*4.5, 'plasmid');
txtl_add_dna(tube3, 'pBAD(50)', 'utr1(20)', 'tetR(600)', vParam(i)*4.5, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'utr1(20)', 'deGFP(1000)-lva(20)',4*4.5, 'plasmid');
txtl_add_dna(tube3,'p70(50)', 'utr1(20)', 'ClpX(600)',0.1*4.5, 'plasmid');

txtl_addspecies(tube3, 'arabinose', 10000);
txtl_addspecies(tube3, 'aTc', 1000);
ClpXToAdd = 80;
% set up well_a1
Mobj{i} = txtl_combine([tube1, tube2, tube3]);

%% Run a simulation

 simulationTime = 0.01*60*60;

simData{i} = txtl_runsim(Mobj{i},simulationTime);

t_ode{i} = simData{i}.Time;
x_ode{i} = simData{i}.Data;

iClpX{i} = findspecies(Mobj{i}, 'protein ClpX');
x_ode{i}(end,iClpX{i}) = x_ode{i}(end,iClpX{i})+ClpXToAdd;%

 simulationTime = 8*60*60;
 [t_ode{i},x_ode{i}] = txtl_runsim(Mobj{i},simulationTime, t_ode{i},x_ode{i});

iGFP{i} = findspecies(Mobj{i}, 'protein deGFP-lva*');
iAraC{i} = findspecies(Mobj{i}, 'protein AraC');
itetR{i} = findspecies(Mobj{i}, 'protein tetRdimer');
iATC{i} = findspecies(Mobj{i}, 'aTc');

end

colororder = lines;

sp1 = zeros(nParam, 1);
sp2 = zeros(nParam, 1);
sp3 = zeros(nParam, 1);
sp4 = zeros(nParam, 1);

figure
subplot(2,2,1)
for i = 1:nParam
sp1(i) = plot(t_ode{i}/60, x_ode{i}(:,iGFP{i}));
set(sp1(i), 'color', colororder(i,:), 'linewidth', 1.5);
hold on
end
legend(sp1, {[Param ' = ' num2str(vParam(1))], [Param ' = ' num2str(vParam(2))], [Param ' = ' num2str(vParam(3))],...
    [Param ' = ' num2str(vParam(4))], [Param ' = ' num2str(vParam(5))]})
title(['GFP conc, ' Param ' varying'])

subplot(2,2,2)
for i = 1:nParam
sp2(i) = plot(t_ode{i}/60, x_ode{i}(:,iAraC{i}));
set(sp2(i), 'color', colororder(i,:), 'linewidth', 1.5);
hold on
end
legend(sp2, {[Param ' = ' num2str(vParam(1)) 'nM'], [Param ' = ' num2str(vParam(2)) 'nM'], [Param ' = ' num2str(vParam(3)) 'nM'],...
    [Param ' = ' num2str(vParam(4)) 'nM'], [Param ' = ' num2str(vParam(5)) 'nM']})
title(['AraC conc, ' Param ' varying'])

subplot(2,2,3)
for i = 1:nParam
sp3(i) = plot(t_ode{i}/60, x_ode{i}(:,itetR{i}));
set(sp3(i), 'color', colororder(i,:), 'linewidth', 1.5);
hold on
end
legend(sp3, {[Param ' = ' num2str(vParam(1))], [Param ' = ' num2str(vParam(2))], [Param ' = ' num2str(vParam(3))],...
    [Param ' = ' num2str(vParam(4))], [Param ' = ' num2str(vParam(5))]})
title(['tetR conc, ' Param ' varying'])

subplot(2,2,4)
for i = 1:nParam
sp4(i) = plot(t_ode{i}/60, x_ode{i}(:,iATC{i}));
set(sp4(i), 'color', colororder(i,:), 'linewidth', 1.5);
hold on
end
legend(sp4, {[Param ' = ' num2str(vParam(1)) 'nM'], [Param ' = ' num2str(vParam(2)) 'nM'], [Param ' = ' num2str(vParam(3)) 'nM'],...
    [Param ' = ' num2str(vParam(4)) 'nM'], [Param ' = ' num2str(vParam(5)) 'nM']})
title(['aTc conc, ' Param ' varying'])

% [ax1,h1]=suplabel('Time/min');
% [ax2,h2]=suplabel('Conc/nM','y');