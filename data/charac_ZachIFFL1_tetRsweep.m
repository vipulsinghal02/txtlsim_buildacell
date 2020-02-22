
close all; clear all


tetRDNAInitial = [0 0.031 0.063 0.125 0.25 0.5 1];
n = length(tetRDNAInitial);
count = 0; Mobj = cell(n,1); simData = cell(n,1); t_ode = cell(n,1); x_ode = cell(n,1); 

for tetRDNA = tetRDNAInitial
    count = count+1;
    tube1 = txtl_extract('E30_1');
    tube2 = txtl_buffer('E30_1');
    tube3 = txtl_newtube('ZachIFFL1');
dna_lasR = txtl_add_dna(tube3,  'placI(50)', 'rbs(20)', 'lasR(1000)',  1, 'plasmid'); % 163 in zach's index	
dna_tetR = txtl_add_dna(tube3,  'plasR(50)', 'SNP10UTR1(20)', 'tetR(1000)', tetRDNA,'plasmid'); % 196
dna_deGFP = txtl_add_dna(tube3,  'plasR_ptet(50)', 'rbs(20)', 'deGFP(1000)-lva', 4,'plasmid');
    
    Mobj{count} = txtl_combine([tube1, tube2, tube3]);
    
txtl_addspecies(Mobj{count}, 'protein ClpX*', 100);
    txtl_addspecies(Mobj{count}, 'OC12HSL', 10000);
    [simData{count}] = txtl_runsim(Mobj{count},12*60*60);
    t_ode{count} = simData{count}.Time;
    x_ode{count} = simData{count}.Data;
%      txtl_plot(simData{count},Mobj{count});
end
figure
colororder = lines;
for i = 1:n
    iGFP = findspecies(Mobj{i}, 'protein deGFP*');
    h(i) = plot(t_ode{i}/60, x_ode{i}(:,iGFP));
    hold on
    set(h(i), 'Color', colororder(i,:), 'LineWidth', 1.5);
    hold on
end        
        title('GFP levels, placI-lasR = 1nM, plasR-tetR = varying, 3OC12HSL = 10uM, plasRptet-GFP = 4nM')  
        legend(h, {'0 nM',' 0.031 nM','0.063 nM','0.125 nM','0.25 nM','0.5 nM','1 nM'}, 'Location', 'NorthEastOutside')
xlabel('time, min')
ylabel('GFP, nM')
        
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
