function mobj = model_protein3(varargin)
% model_protein3 Constitutive gene expression model using a single 
% enzymatic step. 
%
% ~~~ MODEL ~~~
% D + pol <-> D__pol  (k_f, k_r  ) 
% D__pol -> D + pol + protien (kc) 
% 


%% set input defaults
p = inputParser ;
addParameter(p, 'simtime', 1.6*3600);
parse(p);
p = p.Results;

%% setup model
mobj = sbiomodel('expression');

%% setup model reactions
r1 = addreaction(mobj,'dG + pol <-> dG_pol');
Kobj = addkineticlaw(r1,'MassAction');
Kobj.ParameterVariableNames = {'kfdG','krdG'};
addparameter(mobj, 'kfdG', 10);
addparameter(mobj, 'krdG', 600);

r2 = addreaction(mobj,'dG_pol -> dG + pol + pG');
Kobj = addkineticlaw(r2,'MassAction');
Kobj.ParameterVariableNames = {'kcp'};
addparameter(mobj, 'kcp', 0.012);

% setup model species initial concentrations. 
P = sbioselect(mobj, 'name', 'dG');
P.InitialAmount = 0;

C = sbioselect(mobj, 'name', 'pol');
C.InitialAmount = 0;

E = sbioselect(mobj, 'name', 'dG_pol');
E.InitialAmount = 0;

S = sbioselect(mobj, 'name', 'pG');
S.InitialAmount = 0;

%% Run the model

cs = getconfigset(mobj, 'active');
set(cs, 'StopTime', p.simtime);
        
sd = sbiosimulate(mobj);

end

