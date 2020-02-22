% clean up
close all
clc
clearvars

% 
extract = 'E9';

% inducer levels
aTcLevels = [0,5,10,20,50,100,200,500,1000,2000];
IPTGLevels = [0,5,10,20,50,100,200,500,1000,2000];

ptet_DNA = 1;
placI_DNA = 1;

tspan = 10;
IPTG_amount = 0;
aTc_amount = 0;

%% run aTc gradient 
figure(1)
hold on
for k=1:size(aTcLevels,2)
    [simData{k}, Mobj{k}] = genetic_toggle_switch_function(extract,ptet_DNA,placI_DNA,tspan,IPTG_amount,aTcLevels(k));
    tetR = findspecies(Mobj{k},'protein tetRdimer','withInComplex');
    LacI = findspecies(Mobj{k},'protein lacItetramer','withInComplex');
    plot(simData{k}.Time/60,sum(simData{k}.Data(:,tetR),2));
end
hold off

%% IPTG gradient
aTc_amount = 0;

figure(2)
hold on
for k=1:size(IPTGLevels,2)
    [simData{k}, Mobj{k}] = genetic_toggle_switch_function(extract,ptet_DNA,placI_DNA,tspan,IPTGLevels(k),aTc_amount);
    tetR = findspecies(Mobj{k},'protein tetRdimer','withInComplex');
    LacI = findspecies(Mobj{k},'protein lacItetramer','withInComplex');
    plot(simData{k}.Time/60,sum(simData{k}.Data(:,LacI),2));
end
hold off

%% aTc AND IPTG gradient
figure(3)
hold on
for k=1:size(aTcLevels,2)
    for j=1:size(IPTGLevels,2)
        [simData{j}, Mobj{j}] = genetic_toggle_switch_function(extract,ptet_DNA,placI_DNA,tspan,IPTGLevels(j),aTcLevels(k));
        tetR = findspecies(Mobj{j},'protein tetRdimer','withInComplex');
        LacI = findspecies(Mobj{j},'protein lacItetramer','withInComplex');
        plot(aTcLevels(k),simData{j}.Data(end,tetR),'bx');
        plot(IPTGLevels(j),simData{j}.Data(end,LacI),'ro');
    end
    
end
hold off






