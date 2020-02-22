% sbio_genereg.m - simulate and plot output of MATLAB genereg example
% RMM, 7 Sep 2012
%
% This file contains the MATLAB Simbiology gene regulation example, 
% obtains from the online Simbiology documentation.

% Create the TXTL crude extract model
Mobj = sbiomodel('cell');
compObj = addcompartment(Mobj, 'contents');

Robj1 = addreaction(Mobj, '[DNA tetR] -> [DNA tetR] + [mRNA tetR]');
Kobj1 = addkineticlaw(Robj1,'MassAction');
Pobj1 = addparameter(Kobj1, 'k1', 0.2);
set(Kobj1, 'ParameterVariableNames', 'k1');

Sobj1 = sbioselect(Mobj, 'Type', 'species', 'Name', 'DNA tetR');
set(Sobj1, 'InitialAmount', 50);
set(Sobj1, 'UserData', 'ptet(100):rdef(50):tetR(1000)')

Robj2 = addreaction(Mobj, '[mRNA tetR] -> [mRNA tetR] + [protein tetR]');
Kobj2 = addkineticlaw(Robj2,'MassAction');
Pobj2 = addparameter(Kobj2, 'k2', 20.0);
set(Kobj2, 'ParameterVariableNames','k2');

Robj3 = addreaction(Mobj, '[DNA tetR] + [protein tetR] <-> [DNA tetR:protein tetR]');
Kobj3 = addkineticlaw(Robj3,'MassAction');
Pobj3 = addparameter(Kobj3, 'k3', 0.2);
Pobj3r = addparameter(Kobj3, 'k3r', 1.0);
set(Kobj3, 'ParameterVariableNames', {'k3', 'k3r'});

Robj4 = addreaction(Mobj, '[mRNA tetR] -> null');
Kobj4 = addkineticlaw(Robj4, 'MassAction');
Pobj4 = addparameter(Kobj4,  'k4', 1.5);
set(Kobj4, 'ParameterVariableNames','k4');

Robj5 = addreaction(Mobj, '[protein tetR] -> null');
Kobj5 = addkineticlaw(Robj5,'MassAction');
Pobj5 = addparameter(Kobj5,  'k5', 1.0);
set(Kobj5, 'ParameterVariableNames','k5');

% Simulate and plot
[t_ode, x_ode, names] = sbiosimulate(Mobj);
plot(t_ode, x_ode(:,1:4));
legend(names, 'Location', 'NorthEastOutside')
title('Gene Regulation');
xlabel('Time (seconds)');
ylabel('Species Amounts');
