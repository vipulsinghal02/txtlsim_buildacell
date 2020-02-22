% txtl_protein_betaGal.m - protein information for betaGalactosidase
% Zoltan Tuza Sep 2012
%
% This file contains a description of the protein produced by tetR.
% Calling the function txtl_protein_tetR() will set up the reactions for
% sequestration by the inducer aTc.

% Written by Richard Murray, 9 Sep 2012
%
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%
%   2. Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in the 
%      documentation and/or other materials provided with the distribution.
%
%   3. The name of the author may not be used to endorse or promote products 
%      derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
% IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function varargout = txtl_protein_betaGal(mode, tube, protein, varargin)



%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    geneData = varargin{1};
    defaultBasePairs = {'betaGal','lva','terminator';1000,40,100};
    geneData = txtl_setup_default_basepair_length(tube,geneData,...
        defaultBasePairs);
    
    varargout{1} = geneData;
    
    coreSpecies = {'Lactose','alloLactose','Glu+Gal',...
        'protein lacItetramer','alloLactose:protien lacItetramer',...
        'Lactose_ext'};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    
    % Parameters that describe this RBS
% kf_aTc = 1; kr_aTc = 0.1; 

% Set up the binding reaction
%Robj1 = addreaction(tube, [protein.Name ' + aTc <-> aTc:' protein.Name]);
%Kobj1 = addkineticlaw(Robj1, 'MassAction');
%Pobj1f = addparameter(Kobj1, 'kf', kf_aTc);
%Pobj1r = addparameter(Kobj1, 'kr', kr_aTc);
%set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

% 
% % Mass action type version
% % Parameters
% kf_lacbetaGal = 0.6;
% kr_lacbetaGal = 1;
% kf_alloLac = 0.8;
% kf_allolacbetaGalComplex = 0.45;
% kr_allolacbetaGalComplex = 1;
% kf_GluGal = 0.9;
% % conversion of Lactose - enzyme complex
% Robj1 = addreaction(tube, [protein.Name ' + Lactose <-> Lactose:' protein.Name]);
% Kobj1 = addkineticlaw(Robj1, 'MassAction');
% Pobj1f = addparameter(Kobj1, 'kf', kf_lacbetaGal);
% Pobj1r = addparameter(Kobj1, 'kr', kr_lacbetaGal);
% set(Kobj1, 'ParameterVariableNames', {'kf_LacBetaGal', 'kr_LacBetaGal'});
% 
% % forming product
% Robj2 = addreaction(tube, ['Lactose:' protein.Name '-> alloLactose + ' protein.Name]);
% Kobj2 = addkineticlaw(Robj2, 'MassAction');
% Pobj2f = addparameter(Kobj2, 'kf', kf_alloLac);
% set(Kobj2, 'ParameterVariableNames', {'kf_alloLac'});
% 
% % conversion of AlloLactose to Glucose and Galactose
% Robj3 = addreaction(tube, ['alloLactose + ' protein.Name '<-> alloLactose:' protein.Name]);
% Kobj3 = addkineticlaw(Robj3, 'MassAction');
% Pobj3f = addparameter(Kobj3, 'kf', kf_allolacbetaGalComplex);
% Pobj3r = addparameter(Kobj3, 'kr', kr_alolacbetaGalComplex);
% set(Kobj3, 'ParameterVariableNames', {'kf', 'kr'});
% 
% Robj4 = addreaction(tube, ['alloLactose:' protein.Name '-> Glu+Gal+' protein.Name]);
% Kobj4 = addkineticlaw(Robj4, 'MassAction');
% Pobj4f = addparameter(Kobj4, 'kf', kf_GluGal);
% set(Kobj4, 'ParameterVariableNames', {'kf_GluGal'});

    
    % Parameters
    Vl_lac_alloLac = 0.003; %1/min   \alpha_a 
    Kl_lac_alloLac = 1.25;% 
    %! TODO include enzyme concentraction
    Vl_allolac_GluGal = 0.001; %1/min
    Kl_allolac_GluGal =  1.25; %

     kf_seqlacI = 0.5;
     kr_seqlacI = 0.5;

    % conversion of Lactose
    Rlist = addreaction(tube, 'Lactose -> alloLactose');
    Kobj1 = addkineticlaw(Rlist, 'Henri-Michaelis-Menten');
    Pobj1f = addparameter(Kobj1, 'Vl_alloLac',Vl_lac_alloLac);
    Pobj1r = addparameter(Kobj1, 'Kl_alloLac',Kl_lac_alloLac);
    set(Kobj1, 'ParameterVariableNames', {'Vl_alloLac', 'Kl_alloLac'});
    set(Kobj1,'SpeciesVariableNames', {'Lactose'});

    % conversion of AlloLactose to Glucose and Galactose
    Rlist(end+1) = addreaction(tube, 'alloLactose -> Glu+Gal');
    Kobj2 = addkineticlaw(Rlist(end), 'Henri-Michaelis-Menten');
    Pobj2f = addparameter(Kobj2, 'Vl_GluGal',Vl_allolac_GluGal);
    Pobj2r = addparameter(Kobj2, 'Kl_GluGal',Kl_allolac_GluGal);
    set(Kobj2, 'ParameterVariableNames', {'Vl_GluGal', 'Kl_GluGal'});
    set(Kobj2,'SpeciesVariableNames', {'alloLactose'});

    % sequestration of lacI
    % I'm not sure of that it goes here...
    Rlist(end+1) = addreaction(tube,'alloLactose + protein lacItetramer <-> alloLactose:protien lacItetramer');
    Kobj4 = addkineticlaw(Rlist(end), 'MassAction');
    Pobj4f = addparameter(Kobj4, 'kf', kf_seqlacI);
    Pobj4r = addparameter(Kobj4, 'kr', kr_seqlacI);
    set(Kobj4, 'ParameterVariableNames', {'kf','kr'});

    Vl = 0.006;
    Kl = 1.25;
    % transporting Lactose_ext into the cell
    Rlist(end+1)= addreaction(tube, 'Lactose_ext -> Lactose');
    Kobj1 = addkineticlaw(Rlist(end), 'Henri-Michaelis-Menten');
    Pobj1f = addparameter(Kobj1, 'Vl_Lactose_ext',Vl);
    Pobj1r = addparameter(Kobj1, 'Kl_Lactose_ext',Kl);
    set(Kobj1, 'ParameterVariableNames', {'Vl_Lactose_ext', 'Kl_Lactose_ext'});
    set(Kobj1,'SpeciesVariableNames', {'Lactose_ext'});

%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_protein_betaGal:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end    


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
