% protstyn_plot.m - plot results of protein synthesis run
% Clare Chen, Sep 2012

% gamS concentrations for sample run
gamS = [0; 0.02; 0.1; 0.5; 2.5];
szg = size(gamS);
tspan = 0:0.01*3600:1.98*3600;
szt = size(tspan);
graph = zeros(szt(1, 2), 2*szg(1, 1));

for i =1:szg(1,1)
    % initial conditions (initial DNA template concentration of 4 nM)
    init=[4;0;0;0]; 
    chgpar = gamS(i);   
    [t,y]=ode23(@protsynt81bis,tspan,init,[], chgpar);
    graph(:, 2*i-1:2*i) = [t/3600, y(:, 4)/1000];
end

figure(1)
plot(graph(:,1), graph(:,2), '-r',  ...
     graph(:, 3), graph(:, 4), '-b',  ...
     graph(:, 5), graph(:, 6), '-m', ...
     graph(:, 7), graph(:, 8), '-k', ...
     graph(:, 9), graph(:,10), '-c');
hleg = legend('0 然 gamS', '0.02 然 gamS','0.1 然 gamS', ...
              '0.5 然 gamS', '2.5 然 gamS');
title('Simulation of Gene Expression to 2 hours');
xlabel('Time [h]');
ylabel('Scaled Intensity (A.U.)');
