function mobj = model_aTc_induc1
% aTc derepression (induction) of a circuit involving	
% repression with enzymatic one step protein production
% The 1 at the end of the aTc_induc1 is used because this uses
% a transcriptional model of type 1. 
% 
% ~~~ MODEL ~~~
% D_T + P <-> D_T:P -> D_T + P + T
% D_G + P <-> D_G:P -> D_G + P + G
% 2 T <-> T2
% D_G + T2 <-> D_G:T2
% 2 aTc <-> aTc2
% aTc2 + T2 <-> aTc2:T2

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

% D_T + P <-> D_T:P -> D_T + P + T
% D_G + P <-> D_G:P -> D_G + P + G
% 2 T <-> T2
% D_G + T2 <-> D_G:T2
% 2 aTc <-> aTc2
% aTc2 + T2 <-> aTc2:T2
% 

p = inputParser;
addParameter(p, 'simtime', 8*3600);
parse(p);
p = p.Results;
%% setup model
mobj = sbiomodel('aTcInduc1');
%% setup model reactions

rxn = addreaction(mobj,'dT + pol <-> dT_pol');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfdT','krdT'};
addparameter(mobj, 'kfdT', 1);
addparameter(mobj, 'krdT', 60);

rxn = addreaction(mobj,'dT_pol -> dT + pol + pT');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcp'};
addparameter(mobj, 'kcp', 0.12);

rxn = addreaction(mobj,'dG + pol <-> dG_pol');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfdG','krdG'};
addparameter(mobj, 'kfdG', 0.5);
addparameter(mobj, 'krdG', 30);

rxn = addreaction(mobj,'dG_pol -> dG + pol + pG');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcp'};
% already added to model. 

rxn = addreaction(mobj,'2 pT <-> pT2');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfdimTet','krdimTet'};
addparameter(mobj, 'kfdimTet', .2);
addparameter(mobj, 'krdimTet', 4);

rxn = addreaction(mobj,'dG + pT2 <-> dG_pT2');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfseqTet','krseqTet'};
addparameter(mobj, 'kfseqTet', .2);
addparameter(mobj, 'krseqTet', 4);

% rxn = addreaction(mobj,'2 aTc <-> aTc2');
% Kobj = addkineticlaw(rxn,'MassAction');
% Kobj.ParameterVariableNames = {'kfdimaTc','krdimaTc'};
% addparameter(mobj, 'kfdimaTc', .2);
% addparameter(mobj, 'krdimaTc', 4);

rxn = addreaction(mobj,'aTc2 + pT2 <-> aTc2_pT2');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfseqaTc','krseqaTc'};
addparameter(mobj, 'kfseqaTc', .2);
addparameter(mobj, 'krseqaTc', 4);

% setup model species initial concentrations. 
specie = sbioselect(mobj, 'name', 'dT');
specie.InitialAmount = 0.4;

specie = sbioselect(mobj, 'name', 'dG');
specie.InitialAmount = 4;

specie = sbioselect(mobj, 'name', 'pT');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'pG');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'pol');
specie.InitialAmount = 100;

specie = sbioselect(mobj, 'name', 'dT_pol');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'dG_pol');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'pT2');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'dG_pT2');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'aTc2');
specie.InitialAmount = 0;

% specie = sbioselect(mobj, 'name', 'aTc');
% specie.InitialAmount = 1000;

specie = sbioselect(mobj, 'name', 'aTc2_pT2');
specie.InitialAmount = 0;

%% Run the model

cs = getconfigset(mobj, 'active');
set(cs, 'StopTime', p.simtime);

sd = sbiosimulate(mobj);

end

