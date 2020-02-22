% clean up
close all
clc
clearvars

% 
extract = 'E9';
tspan = 7; % hours

% inducer levels
aTcLevels = [0,5,10,20,50,100,200,500,1000,2000];
arabinoseLevels = [0,5,10,20,50,100,200,500,1000,2000];
colors = {'r', 'b', 'g', 'c', 'm', 'y', 'k', 'r--', 'b--','c--'};


p70_AraC = 1;
pBAD_tetR = 1.1;
pBAD_ptet_deGFP = 1.1;
arabinose_amount = 1000;
aTc_amount = 1;
protein_ClpX =100;

%% run aTc gradient 
figure(1)
hold on
for k=1:size(aTcLevels,2)
    [simData{k}, Mobj{k}] = incoherent_ff_loop_function(extract,p70_AraC,pBAD_tetR,pBAD_ptet_deGFP,tspan,protein_ClpX,arabinose_amount,aTcLevels(k));
    
    GFP = findspecies(Mobj{k},'protein deGFP-lva*');
    plot(simData{k}.Time/60,sum(simData{k}.Data(:,GFP),2),colors{k});
    labels{k} = [int2str(aTcLevels(k)) ' nM aTc'];
end
xlabel('Time [min]')
ylabel('GFP-lva* [nM]')
legend(labels)
hold off



%% run aTc gradient 
figure(2)
hold on
for k=1:size(aTcLevels,2)
    [simData{k}, Mobj{k}] = incoherent_ff_loop_function(extract,p70_AraC,pBAD_tetR,pBAD_ptet_deGFP,tspan,protein_ClpX,arabinoseLevels(k),aTc_amount);
    
    GFP = findspecies(Mobj{k},'protein deGFP-lva*');
    plot(simData{k}.Time/60,sum(simData{k}.Data(:,GFP),2),colors{k});
    labels{k} = [int2str(aTcLevels(k)) ' nM arabinose'];
end
xlabel('Time [min]')
ylabel('GFP-lva* [nM]')
legend(labels)
hold off

