%% clean up

clear variables
clc
close all
%% My test run

% Set up the standard TXTL tubes
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('E6');

% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
dna_lacI = txtl_adddna(tube3, ...
    'p70(50)', 'rbs(20)', 'lacI(600)', 3, 'linear');
dna_deGFP = txtl_adddna(tube3, ...
    'placI(50)', 'rbs(20)', 'deGFP(1000)', 3, 'linear');
dna_gamS = txtl_adddna(tube3, ...
    'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes and add some inducer
well_a1 = txtl_combine([tube1, tube2, tube3]);


%% Run a simulation

simulationTime = 4.5*60*60;


% 1st run
[t_ode,x_ode] = txtl_runsim(well_a1,simulationTime);

disp('1')
pause(1)
% 2nd run
%txtl_continue_simulation(simData,well_a1);

% add more DNA
txtl_addspecies(well_a1, 'DNA p70--rbs--lacI', 1);
txtl_addspecies(well_a1, 'DNA placI--rbs--deGFP', 1);

[t_ode_2,x_ode_2] = txtl_runsim(well_a1,simulationTime,t_ode, x_ode);
disp('2')
pause(1)


% 3rd run
%txtl_continue_simulation(simData_2,well_a1);

% add more DNA
txtl_addspecies(well_a1, 'DNA p70--rbs--lacI', 1);
txtl_addspecies(well_a1, 'DNA placI--rbs--deGFP', 1);


[t_ode_3,x_ode_3] = txtl_runsim(well_a1,simulationTime,t_ode_2,x_ode_2);
disp('3')
pause(1)
% concatante data


%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'};
dataGroups{1,3} = {'b-','r-','b--','r--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*','protein gamS','[protein lacI]_tot','[protein gamS]_tot'};
dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,well_a1,dataGroups)
