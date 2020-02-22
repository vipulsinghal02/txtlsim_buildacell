
close all; clear all


inducerInitial = [0.01 0.1 1 10 100 1000];
n = length(inducerInitial);
count = 0; Mobj = cell(n,1); simData = cell(n,1); t_ode = cell(n,1); x_ode = cell(n,1); 

for inducer = inducerInitial
    count = count+1;
    tube1 = txtl_extract('E30_1');
    tube2 = txtl_buffer('E30_1');
    tube3 = txtl_newtube('gene_expression');

    dna_tetR = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'lasR(1000)',  1,  'plasmid');					% type
    dna_deGFP = txtl_add_dna(tube3, 'plas(50)', 'rbs(20)', 'deGFP(1000)', 2, 'plasmid');
    
    Mobj{count} = txtl_combine([tube1, tube2, tube3]);
    txtl_addspecies(Mobj{count}, 'OC12HSL', inducer);
    [simData{count}] = txtl_runsim(Mobj{count},8*60*60);
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
        title('GFP levels, p70-lasR = 1nM, plasR-GFP = 2nM, 3OC12HSL varying')  
        legend(h, {'0.01 nM', '0.1 nM', '1 nM', '10 nM', '100 nM', '1000 nM'}, 'Location', 'NorthEastOutside')
    figure
    finalGFP = zeros(1,n);
for i = 1:n
    iGFP = findspecies(Mobj{i}, 'protein deGFP*');
    finalGFP(i) = x_ode{i}(end,iGFP);
end        

        semilogx(inducerInitial', finalGFP', '*-', 'LineWidth', 2)
        title('GFP endpoint values')  
        xlabel('OC12HSL, nM')
        ylabel('GFP endpoint, nM')
        
        
        
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
