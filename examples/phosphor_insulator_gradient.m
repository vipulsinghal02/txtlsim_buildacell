% clean up
close all
clc
clearvars

% 
extract = 'E30VNPRL';
tspan = 7; % hours

% inducer levels
NRIDNALeves = [0,5,10,20,50,100,200,500,1000,2000];
colors = {'r', 'b', 'g', 'c', 'm', 'y', 'k', 'r--', 'b--','c--'};


dna_sigma54 = 1;
dna_NRII_L16R = 1;
dna_NRII_H139N = 1;
dna_deGFP = 1;


%% run aTc gradient 
figure(1)
hold on
for k=1:size(NRIDNALeves,2)
    [simData{k}, Mobj{k}] = phosphor_insulator_function(extract,dna_sigma54,dna_NRII_L16R,dna_NRII_H139N,NRIDNALeves(k),dna_deGFP,tspan);
    
    NRI = findspecies(Mobj{k},'protein NRI','withInComplex');
    plot(simData{k}.Time/60,sum(simData{k}.Data(:,NRI),2),colors{k});
    labels{k} = [int2str(NRIDNALeves(k)) ' nM NRI dna'];
end
xlabel('Time [min]')
ylabel('NRI [nM]')
legend(labels)
hold off





