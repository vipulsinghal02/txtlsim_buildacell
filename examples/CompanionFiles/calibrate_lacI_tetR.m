% Calibrate lacI and tetR components
close all
clc
clear all
B_E = 'E9';
SimulationNumber = '9';
mainDir = pwd;

%% COnstitutive Production
tube1 = txtl_extract(B_E);
tube2 = txtl_buffer(B_E);
tube3 = txtl_newtube('GE_ptet');
tube4 = txtl_newtube('GE_ptet3');
tube5 = txtl_newtube('GE_placI');
tube6 = txtl_newtube('GE_placI3');
dna_deGFP = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'tetR(1000)',1,'plasmid');
dna_deGFP = txtl_add_dna(tube4, 'p70(50)', 'rbs(20)', 'tetR3(1000)',1,'plasmid');
dna_deGFP = txtl_add_dna(tube5, 'p70(50)', 'rbs(20)', 'lacI(1000)',1,'plasmid');
dna_deGFP = txtl_add_dna(tube6, 'p70(50)', 'rbs(20)', 'lacI3(1000)',1,'plasmid');
Mobj1 = txtl_combine([tube1, tube2, tube3]);
Mobj2 = txtl_combine([tube1, tube2, tube4]);
Mobj3 = txtl_combine([tube1, tube2, tube5]);
Mobj4 = txtl_combine([tube1, tube2, tube6]);

figure
itetR = findspecies(Mobj1, 'protein tetR');
itetRdimer = findspecies(Mobj1, 'protein tetRdimer');
[simData1] = txtl_runsim(Mobj1,14*60*60);
t_ode1 = simData1.Time;
x_ode1 = simData1.Data;
p1 = plot(t_ode1, x_ode1(:,itetR), 'b', t_ode1, x_ode1(:,itetRdimer), 'r')
hold on

[simData2] = txtl_runsim(Mobj2,14*60*60);
t_ode2 = simData2.Time;
x_ode2 = simData2.Data;
itetR = findspecies(Mobj2, 'protein tetR3');
itetRdimer = findspecies(Mobj2, 'protein tetR3dimer');
p2 = plot(t_ode2, x_ode2(:,itetR), 'b--', t_ode2, x_ode2(:,itetRdimer), 'r--')
legend([p1;p2], 'tetR', 'tetRdimer', 'tetR3', 'tetR3dimer', 'Location', 'NorthEastOutside')
title('tetR constitutive Production')
xlabel('time/h'); ylabel('conc/nM');
cd([mainDir '\examples\Characterization work\Characterizing ptet tetR placI lacI'])
print('-dtiff','-r80',['Sim_' SimulationNumber '_1_'])
cd(mainDir)


figure
[simData3] = txtl_runsim(Mobj3,14*60*60);
t_ode3 = simData3.Time;
x_ode3 = simData3.Data;
ilacI = findspecies(Mobj3, 'protein lacI');
ilacIdimer = findspecies(Mobj3, 'protein lacIdimer');
ilacItetramer = findspecies(Mobj3, 'protein lacItetramer');
p3 = plot(t_ode3, x_ode3(:,ilacI), 'b', t_ode3, x_ode3(:,ilacIdimer), 'r',  t_ode3, x_ode3(:,ilacItetramer), 'g')
hold on
[simData4] = txtl_runsim(Mobj4,14*60*60);
t_ode4 = simData4.Time;
x_ode4 = simData4.Data;
ilacI = findspecies(Mobj4, 'protein lacI3');
ilacIdimer = findspecies(Mobj4, 'protein lacI3dimer');
ilacItetramer = findspecies(Mobj4, 'protein lacI3tetramer');
p4 = plot(t_ode4, x_ode4(:,ilacI), 'c--', t_ode4, x_ode4(:,ilacIdimer), 'm--',  t_ode4, x_ode4(:,ilacItetramer), 'k--')

legend([p3;p4], 'lacI', 'lacIdimer', 'lacItetramer', 'lacI3', 'lacI3dimer',...
    'lacI3tetramer', 'Location', 'NorthEastOutside')
title('lacI constitutive Production')
xlabel('time/h'); ylabel('conc/nM');
cd([mainDir '\examples\Characterization work\Characterizing ptet tetR placI lacI'])
print('-dtiff','-r80',['Sim_' SimulationNumber '_2_'])
cd(mainDir)

%% Constitutive Production Type 2
tube1 = txtl_extract(B_E);
tube2 = txtl_buffer(B_E);
tube3 = txtl_newtube('GE2_ptet3');
tube4 = txtl_newtube('GE2_placI3');

dna_GFP = txtl_add_dna(tube3, 'ptet3(50)', 'rbs(20)', 'deGFP(1000)',1,'plasmid');
dna_RFP = txtl_add_dna(tube4, 'placI3(50)', 'rbs(20)', 'RFP(1000)',1,'plasmid');
Mobj1 = txtl_combine([tube1, tube2, tube3]);
Mobj2 = txtl_combine([tube1, tube2, tube4]);

figure
% Dont need the FP from ptet because ptet and ptet3 are exactly the same in
% their TX and TL rates. Theu differ in DNA sequestration and multimerization
% of the corresponding protein), which is not present in Constitutive production of FP)

[simData1] = txtl_runsim(Mobj1,14*60*60);
t_ode1 = simData1.Time;
x_ode1 = simData1.Data;
iGFP = findspecies(Mobj1, 'protein deGFP*');
p1 = plot(t_ode1, x_ode1(:,iGFP), 'g')
hold on

[simData2] = txtl_runsim(Mobj2,14*60*60);
t_ode2 = simData2.Time;
x_ode2 = simData2.Data;
iRFP= findspecies(Mobj2, 'protein RFP*');
p2 = plot(t_ode2, x_ode2(:,iRFP), 'r--')

legend([p1;p2], 'GFP* (ptet3)', 'RFP* (placI3)', 'Location', 'NorthEastOutside')
title('constitutive Production, ptet3-GFP, placI3-RFP')
xlabel('time/h'); ylabel('conc/nM');
cd([mainDir '\examples\Characterization work\Characterizing ptet tetR placI lacI'])
print('-dtiff','-r80',['Sim_' SimulationNumber '_5_'])
cd(mainDir)

%% Negative Autoregulation
tube3 = txtl_newtube('NA_ptet');
tube4 = txtl_newtube('NA_ptet3');
tube5 = txtl_newtube('NA_placI');
tube6 = txtl_newtube('NA_placI3');

dna_tetR = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'tetR(1000)',1,'plasmid');
dna_deGFP = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'deGFP(1000)',1,'plasmid');
dna_tetR3 = txtl_add_dna(tube4, 'ptet3(50)', 'rbs(20)', 'tetR3(1000)',1,'plasmid');
dna_CFP = txtl_add_dna(tube4, 'ptet3(50)', 'rbs(20)', 'RFP(1000)',1,'plasmid');
dna_lacI = txtl_add_dna(tube5, 'placI(50)', 'rbs(20)', 'lacI(1000)',1,'plasmid');
dna_deGFP = txtl_add_dna(tube5, 'placI(50)', 'rbs(20)', 'deGFP(1000)',1,'plasmid');
dna_lacI3 = txtl_add_dna(tube6, 'placI3(50)', 'rbs(20)', 'lacI3(1000)',1,'plasmid');
dna_CFP = txtl_add_dna(tube6, 'placI3(50)', 'rbs(20)', 'RFP(1000)',1,'plasmid');

Mobj1 = txtl_combine([tube1, tube2, tube3]);
Mobj2 = txtl_combine([tube1, tube2, tube4]);
Mobj3 = txtl_combine([tube1, tube2, tube5]);
Mobj4 = txtl_combine([tube1, tube2, tube6]);

figure
itetR = findspecies(Mobj1, 'protein tetR');
itetRdimer = findspecies(Mobj1, 'protein tetRdimer');
iGFP = findspecies(Mobj1, 'protein deGFP*');
[simData1] = txtl_runsim(Mobj1,14*60*60);
t_ode1 = simData1.Time;
x_ode1 = simData1.Data;
p1 = plot(t_ode1, x_ode1(:,itetR), 'b', t_ode1, x_ode1(:,itetRdimer), 'k', t_ode1, x_ode1(:,iGFP), 'g')
hold on
[simData2] = txtl_runsim(Mobj2,14*60*60);
t_ode2 = simData2.Time;
x_ode2 = simData2.Data;
itetR = findspecies(Mobj2, 'protein tetR3');
itetRdimer = findspecies(Mobj2, 'protein tetR3dimer');
iRFP = findspecies(Mobj2, 'protein RFP*');
p2 = plot(t_ode2, x_ode2(:,itetR), 'b--', t_ode2, x_ode2(:,itetRdimer), 'k--',t_ode2, x_ode2(:,iRFP), 'r')
legend([p1;p2], 'tetR', 'tetRdimer','deGFP*', 'tetR3', 'tetR3dimer', 'RFP*',...
    'Location', 'NorthEastOutside')
title('tetR negative autoregulation')
xlabel('time/h'); ylabel('conc/nM');
cd([mainDir '\examples\Characterization work\Characterizing ptet tetR placI lacI'])
print('-dtiff','-r80',['Sim_' SimulationNumber '_3_'])
cd(mainDir)

figure
[simData3] = txtl_runsim(Mobj3,14*60*60);
t_ode3 = simData3.Time;
x_ode3 = simData3.Data;
ilacI = findspecies(Mobj3, 'protein lacI');
ilacIdimer = findspecies(Mobj3, 'protein lacIdimer');
ilacItetramer = findspecies(Mobj3, 'protein lacItetramer');
iGFP = findspecies(Mobj3, 'protein deGFP*');
p3 = plot(t_ode3, x_ode3(:,ilacI), 'b', t_ode3, x_ode3(:,ilacIdimer), 'b--',...
    t_ode3, x_ode3(:,ilacItetramer), 'b-.',  t_ode3, x_ode3(:,iGFP), 'g')
hold on

[simData4] = txtl_runsim(Mobj4,14*60*60);
t_ode4 = simData4.Time;
x_ode4 = simData4.Data;
ilacI = findspecies(Mobj4, 'protein lacI3');
ilacIdimer = findspecies(Mobj4, 'protein lacI3dimer');
ilacItetramer = findspecies(Mobj4, 'protein lacI3tetramer');
iRFP = findspecies(Mobj4, 'protein RFP*');
p4 = plot(t_ode4, x_ode4(:,ilacI), 'k', t_ode4, x_ode4(:,ilacIdimer), 'k--',...
    t_ode4, x_ode4(:,ilacItetramer), 'k-.',t_ode4, x_ode4(:,iRFP), 'r')
legend([p3;p4], 'lacI', 'lacIdimer', 'lacItetramer','deGFP*', 'lacI3', ...
    'lacI3dimer', 'lacI3tetramer','RFP*', 'Location', 'NorthEastOutside')
title('lacI negative autoregulation')
xlabel('time/h'); ylabel('conc/nM');
cd([mainDir '\examples\Characterization work\Characterizing ptet tetR placI lacI'])
print('-dtiff','-r80',['Sim_' SimulationNumber '_4_'])
cd(mainDir)







