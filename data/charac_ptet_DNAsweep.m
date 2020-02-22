
close all; clear all



count = 0; Mobj = cell(5,1); simData = cell(5,1); t_ode = cell(5,1); x_ode = cell(5,1); 
for inducerInitial = 1000*[0.001 0.01 0.1 1 10]
    count = count+1;
    tube1 = txtl_extract('E30_1');
    tube2 = txtl_buffer('E30_1');
    tube3 = txtl_newtube('gene_expression');

    dna_tetR = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'tetR(1000)',  1,  'plasmid');					% type
    dna_deGFP = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'deGFP(1000)', 2, 'plasmid');
    
    Mobj{count} = txtl_combine([tube1, tube2, tube3]);
    txtl_addspecies(Mobj{count}, 'aTc', inducerInitial);
    [simData{count}] = txtl_runsim(Mobj{count},8*60*60);
    t_ode{count} = simData{count}.Time;
    x_ode{count} = simData{count}.Data;
%      txtl_plot(simData{count},Mobj{count});
end
figure
colororder = lines;
for i = 1:5
    iGFP = findspecies(Mobj{i}, 'protein deGFP*')
    h(i) = plot(t_ode{i}/60, x_ode{i}(:,iGFP));
    hold on
    set(h(i), 'Color', colororder(i,:), 'LineWidth', 1.5);
    hold on
end        
        title('GFP levels, p70-tetR = 1nM, ptet-GFP = 2nM, aTc varying')  
        legend(h, {'1', '10', '100', '1000', '10000'}, 'Location', 'NorthEastOutside')
    figure
    finalGFP = zeros(1,5);
for i = 1:5
    iGFP = findspecies(Mobj{i}, 'protein deGFP*');
    finalGFP(i) = x_ode{i}(end,iGFP);
end        
        semilogx(1000*[0.001 0.01 0.1 1 10], finalGFP)
        title('GFP endpoint values')  
        xlabel('aTc, nM')
        ylabel('GFP endpoint, nM')
        
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
