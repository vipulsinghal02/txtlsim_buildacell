function mobj = model_protein5
% enzymatic mrna and protein production and first order mrna degradation
% 
% ~~~ MODEL ~~~
% D + pol <-> D__pol  (k_fd, k_rd) 
% D__pol -> D + pol + mrna (kcm) 
% 
% mrna + ribo <-> mrna__ribo (k_fm, k_rm)
% mrna__ribo <-> mrna + ribo + protein (kcp)
% 
% mrna -> null (kcx)
% 
% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

p = inputParser ;
addParameter(p, 'simtime', 2*3600);
parse(p);
p = p.Results;


% Model Parameter values
cpol = 100; % nM
cribo = 50; %nM

rkfdG = 10; % nM-1s-1
rkrdG = 600; % s-1
rkcm = 0.001; %s-1

rkfpG = 10; % nM-1s-1
rkrpG = 300; % s-1
rkcp = 1/36;

rdel_m = log(2)/720; % 12 min half life of mrna

%% setup model reactions

% setup model
mobj = sbiomodel('expression');
% GFP TXTL
rxn = addreaction(mobj,'dG + pol <-> dG_pol');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfdG','krdG'};
addparameter(mobj, 'kfdG', rkfdG);
addparameter(mobj, 'krdG', rkrdG);

rxn = addreaction(mobj,'dG_pol -> dG + pol + mG');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcm'};
addparameter(mobj, 'kcm', rkcm);

rxn = addreaction(mobj,'mG + ribo <-> mG_ribo');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfpG','krpG'};
addparameter(mobj, 'kfpG', rkfpG);
addparameter(mobj, 'krpG', rkrpG);

rxn = addreaction(mobj,'mG_ribo -> mG + ribo + pG');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcp'};
addparameter(mobj, 'kcp', rkcp);

rxn = addreaction(mobj,'mG -> null');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'del_m'};
addparameter(mobj, 'del_m', rdel_m);

% setup model species initial concentrations. 

specie = sbioselect(mobj, 'name', 'dG');
specie.InitialAmount = 30;

specie = sbioselect(mobj, 'name', 'pG');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'mG');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'pol');
specie.InitialAmount = cpol;

specie = sbioselect(mobj, 'name', 'ribo');
specie.InitialAmount = cribo;

specie = sbioselect(mobj, 'name', 'dG_pol');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'mG_ribo');
specie.InitialAmount = 0;


%% Run the model

cs = getconfigset(mobj, 'active');
set(cs, 'StopTime', p.simtime);
        
sd = sbiosimulate(mobj);

end

