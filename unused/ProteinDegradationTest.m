% Test for Protein Degradation

tube = txtl_newtube('proteintube');
protein = addspecies(tube, 'protein tetR-lva-terminator', 30);
protein.UserData = 500 / 3;
degradationRate = [0.001 0.01 0.1]; % !TODO: Find a reasonable value
Rlist = txtl_protein_degradation(tube, protein,degradationRate);

configsetObj = getconfigset(tube, 'active');
set(configsetObj, 'StopTime', 5*60*60)
if ~strcmp(version('-release'),'2012a')
 set(configsetObj, 'SolverType', 'ode23s');
end
[t_ode, x_ode, names] = sbiosimulate(tube, configsetObj);

figure(1); clf();
iTetR = findspecies(tube, 'protein tetR-lva-terminator')
p = plot(t_ode/60, x_ode(:, iTetR), 'b-');
set(p, 'LineWidth',2);
set(gca, 'FontSize', 14);
title('Protein Degradation Example');
lgh = legend({'TetR'}, 'Location', 'NortheastOutside');
%legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
axis([0 100 0 30]);
tube.Species