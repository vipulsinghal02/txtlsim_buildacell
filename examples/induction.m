% induction.m - negative autoregulation with an inducer
% R. M. Murray, 11 Sep 2012

%    http://www.mathworks.com/help/toolbox/simbio/gs/fp58748.html
%
clear vars
close all
clear all
%% Do runs at different inducer levels, linearly spaced
levels = [0 2 5 10 20 40 100 300 1000];
maxtetR = zeros(1, length(levels));
colors = {'r', 'b', 'g', 'c', 'm', 'y', 'k', 'r--', 'b--'};

tspan = 10; % 10 hours
atc_amount = [0 2 5 10 20 40 100 300 1000];
dna_amount = 1; % in nM

figure(1)
hold on
for k = 1:size(levels,2)
    [t_ode, x_ode,Mobj] = negautoreg_function('E9',dna_amount,tspan,atc_amount(k));
    tetR = findspecies(Mobj,'protein tetRdimer','withInComplex');
    plot(t_ode/60,sum(x_ode(:,tetR),2),colors{k})
    
    % Keep track of the max expression for later plotting
    maxtetR(k) = sum(x_ode(end,tetR),2);
    labels{k} = [int2str(atc_amount(k)) ' nM aTc'];
end
hold off

% Label the time trace
title('Time Responses');
lgh = legend(labels, 'Location', 'NorthWest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

% Plot the characteristic curve
figure(2); clf();
title('Characteristic Curve');
plot(levels, maxtetR);
xlabel('aTc concentration [nM]');
ylabel('max tetR expression [nM]');



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
