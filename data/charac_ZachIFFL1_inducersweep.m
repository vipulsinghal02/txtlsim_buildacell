
close all; clear all


inducerInitial = [0 0.1 1 10 100 1000 10000];
n = length(inducerInitial);
count = 0; Mobj = cell(n,1); simData = cell(n,1); t_ode = cell(n,1); x_ode = cell(n,1); 

for inducer = inducerInitial
    count = count+1;
    tube1 = txtl_extract('E30_1');
    tube2 = txtl_buffer('E30_1');
    tube3 = txtl_newtube('ZachIFFL1');
dna_lasR = txtl_add_dna(tube3,  'placI(50)', 'rbs(20)', 'lasR(1000)',  1, 'plasmid'); % 163 in zach's index	
dna_tetR = txtl_add_dna(tube3,  'plasR(50)', 'SNP10UTR1(20)', 'tetR(1000)', 0.005,'plasmid'); % 196
dna_deGFP = txtl_add_dna(tube3,  'plasR_ptet(50)', 'rbs(20)', 'deGFP(1000)', 4,'plasmid');
    
    Mobj{count} = txtl_combine([tube1, tube2, tube3]);
    txtl_addspecies(Mobj{count}, 'OC12HSL', inducer);
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
        title('GFP levels, placI-lasR = 1nM, plasR-tetR = 0.005 nM, 3OC12HSL = varying, plasRptet-GFP = 4nM')  
        legend(h, {'0 nM',' 0.1 nM','1 nM','10 nM','100 nM','1000 nM','10000 nM'}, 'Location', 'NorthEastOutside')
xlabel('time, min')
ylabel('GFP, nM')
        
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
