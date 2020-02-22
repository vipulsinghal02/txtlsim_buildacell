%TXTL_ADD_DNA   Set up species and reactions for a DNA segment
%
%
%   dna = TXTL_ADD_DNA(tube, prom_spec, utr_spec, cds_spec, amount, type)
%   constructs the species and reactions required for transcription,
%   translation and degradation of DNA, mRNA and proteins in the
%   TX-TL system.
%
%   * tube = Simbiology model object
%   * preprom_spec = Cell array of nucleatide sequences and corresponding
%   sizes. One example of their use is as a protection from exonucleases.
%   * prom_spec = spec of the form 'pre_prom(nn)'-'prom(nn)' where 'prom' is the
%     promoter name and 'len' is the length of the promoter. pre_prom cound
%     consist of nucleatide sequences and corresponding
%   sizes. One example of their use is as a protection from exonucleases.
%   * utr_spec = spec of the form 'rbs(nn)' where 'rbs' is the RBS
%     name and 'len' is the length of the RBS.
%   * cds_spec = spec of the form 'gene(nn)-lva(nn)-terminator(nn)' where 'gene' is the
%     gene name and 'len' is the length of the gene.
%   * amount = amount of DNA to put in the tube (in nM)
%   * type = 'linear' if you want to include degradation reactions

% Written by Richard Murray, Sep 2012
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
%%
function dna = txtl_add_dna(tube, prom_spec, utr_spec, cds_spec, dna_amount, type, varargin)
mode = struct('add_dna_driver', {[]},...
              'prot_deg_flag',{false},...
              'no_protein_flag',{false}, ...
              'prot_term_flag',{false},...
              'prom_junk_flag',{false},...
              'prom_thio_flag',{false}); 
          
% Extract out the names and lengths
[promData, promStr] = txtl_parsespec(prom_spec);
[utrData, utrStr] = txtl_parsespec(utr_spec);
[geneData, geneStr] = txtl_parsespec(cds_spec);

% check for variations in DNA, used to select specific code
mode.prot_deg_flag = checkForStringInACellList(geneData(1,:),'lva');
mode.prot_term_flag = checkForStringInACellList(geneData(1,:),'terminator');
mode.no_protein_flag = checkForStringInACellList(geneData(1,:),'no_protein');
[mode.prom_junk_flag, junkIndex] = checkForStringInACellList(promData(1,:),'junk');
mode.prom_thio_flag = checkForStringInACellList(promData(1,:),'thio');

% species name string building
geneName = geneData{1,1}; %assuming the format is gene-lva-...-terminator

rbsName = utrData{1,end}; % format is att-...-rbs.
if mode.no_protein_flag
    protstr = ['protein ' geneStr];
    rnastr = ['RNA ' utrStr];
    dnastr = ['DNA ' promStr '--' utrStr];
else
    protstr = ['protein ' geneStr]; % protstr looks something like 'protein tetR-lva-terminator'
    rnastr = ['RNA ' utrStr '--' geneStr];
    dnastr = ['DNA ' promStr '--' utrStr '--' geneStr];
end

promoterName = promData{1,end}; % assuming {'thio','junk','prom'}




%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    mode.add_dna_driver = 'Setup Species';
    tubeUser = get(tube, 'UserData');
    dnaList = tubeUser.DNAinfo;
    dnaList{end+1} = {prom_spec, utr_spec, cds_spec, dna_amount, type, 'rxns_not_set_up', mode};
    tubeUser.DNAinfo = dnaList;
    set(tube, 'UserData', tubeUser)
    clear dnaList tubeUser
    
    % promoter reactions require the lengths of the rna, therefore need to be
    % set up after the protein and utr files are called.
    
    %% Protein properties, parameters and reactions %
    protein = txtl_addspecies(tube, protstr, 0, 'Internal');
    if exist(['txtl_protein_' geneName], 'file') == 2
        geneData = eval(['txtl_protein_' geneName '(mode, tube, protein, geneData)']);
    end
    genelenTot = sum(cell2mat(geneData(2,:)));
    protein.UserData = genelenTot / 3;
    
    %% Untranslated Region properties, parameters and reactions %
    
    rna = txtl_addspecies(tube, rnastr, 0, 'Internal');
    if exist(['txtl_5UTR_' rbsName], 'file') == 2
            %ribosome must bind
            [Ribobound, utrlen] = eval(['txtl_5UTR_' rbsName '(mode, tube, rna, protein, utrData)']);
    else 
            warning('txtltoolbox:txtl_add_dna:fileNotFound', ['TXTL: can''t find txtl_5UTR_' rbsName ...
                '; using default rbs params']);
            [Ribobound, utrlen] = txtl_5UTR_utr1(mode, tube, rna, protein, utrData);
    end
    
    utrlenTot = sum(cell2mat(utrlen(2,:)));
    RNAlengthData = utrlenTot + genelenTot;
    rna.UserData = RNAlengthData;
    
    %% Promoter properties, parameters and reactions %%%%%%%%%%%%%%%%%%%%%%
    % DNA solution is 22.5% of the 10ul reaction volume
    stockMulti = 10/2.25;
    dna_amount = dna_amount*stockMulti;
    dna = txtl_addspecies(tube, dnastr, dna_amount, 'Internal');
    % Transcription %%
    if exist(['txtl_prom_' promoterName], 'file') == 2 
        promData = eval(['txtl_prom_' promoterName '(mode, tube, dna, rna, promData, prom_spec, utr_spec, cds_spec)']);
    else
        warning(['TXTL: can''t find txtl_prom_' promoterName ...
            '; using default promoter params']);
        promData = txtl_prom_p70(mode, tube, dna, rna, promData);
    end
    promlenTot = sum(cell2mat(promData(2,:)));
    %     if mode.utr_antisense_flag
    %         dna.UserData = promlenTot + utrlenTot;
    %     else
    %         dna.UserData = promlenTot + utrlenTot + genelenTot;
    %     end
    dna.UserData = promlenTot + utrlenTot + genelenTot;
    
    %% Translation %%
         txtl_translation(mode, tube, dna, rna, protein, Ribobound);
    
    
    %% DNA, protein degradation
    
    % DNA degradation
    if strcmp(type, 'linear') 
        txtl_dna_degradation(mode, tube, dna);
    end
    
    % Protein degradation (if tagged)
    if mode.prot_deg_flag
        txtl_protein_degradation(mode, tube, protein);
    end
    
    % All done!
    return
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(varargin{1}, 'Setup Reactions')
    mode.add_dna_driver = varargin{1};
    % get a list of the species to search through before setting up
    % certain reactions
    listOfSpecies = get(tube.species, 'name');
    
    
    % set up protein reactions and data, followed by utr followed by promoter
    % (promoter reactions require the lengths of the rna, therefore need to be
    % set up after the protein and utr files are called.
    
    %% Protein properties, parameters and reactions %%%%%%%%%%%%%%%%%%%%%%%
    
    protein = sbioselect(tube, 'Name', protstr);
    
    if exist(['txtl_protein_' geneName], 'file') == 2
        % Run the protein specific setup
        eval(['txtl_protein_' geneName '(mode, tube, protein, listOfSpecies)']);
    end
    
    %% Untranslated Region properties, parameters and reactions %%%%%%%%%%%
    
    rna = sbioselect(tube, 'Name', rnastr);
    if exist(['txtl_5UTR_' rbsName], 'file') == 2
        eval(['txtl_5UTR_' rbsName '(mode, tube, rna, protein)']);
    else
        txtl_5UTR_utr1(mode, tube, rna, protein, utrData);
    end
    
    %% Promoter properties, parameters and reactions %%%%%%%%%%%%%%%%%%%%%%
    
    dna = sbioselect(tube, 'Name', dnastr);
    % Transcription %%
    if exist(['txtl_prom_' promoterName], 'file') == 2 
        eval(['txtl_prom_' promoterName '(mode, tube, dna, rna, listOfSpecies,prom_spec, utr_spec, cds_spec)']);

    else
        warning(['TXTL: can''t find txtl_prom_' promoterName ...
            '; using default promoter params']);
        txtl_prom_p70(mode, tube, dna, rna, listOfSpecies);
    end
    
    % Translation %%
        Ribobound = sbioselect(tube, 'Name', ['Ribo:' rna.Name]);
        txtl_translation(mode, tube, dna, rna, protein, Ribobound);
    %% DNA, mRNA, protein degradation
    
    % DNA degradation
    if strcmp(type, 'linear')
        if mode.prom_junk_flag
            junkLength = promData{2,junkIndex};
            kDNA_complex_deg = log(2)/(1+junkLength/100);
        else
            kDNA_complex_deg = tube.UserData.ReactionConfig.DNA_RecBCD_complex_deg;
        end
        if mode.prom_thio_flag
            kDNA_complex_deg = 0.5*kDNA_complex_deg;
        end
        
        % forward rr for DNA + RecBCD <-> DNA:RecBCD
        kDNA_recbcd_f = tube.UserData.ReactionConfig.DNA_RecBCD_Forward;
        % backward rr for DNA + RecBCD <-> DNA:RecBCD
        kDNA_recbcd_r = tube.UserData.ReactionConfig.DNA_RecBCD_Reverse;
        
        % linear DNA protection by gamS
        txtl_addreaction(tube,'RecBCD + [protein gamS] <-> RecBCD:gamS',...
            'MassAction',{'GamS_RecBCD_f',tube.UserData.ReactionConfig.GamS_RecBCD_F;...
            'GamS_RecBCD_r',tube.UserData.ReactionConfig.GamS_RecBCD_R});
        
        txtl_dna_degradation(mode, tube, dna, [kDNA_recbcd_f, kDNA_recbcd_r, kDNA_complex_deg]);
    end
    
    % Add in mRNA degradation reactions
    txtl_mrna_degradation(mode, tube, dna, rna, utr_spec);
    
    % Protein degradation (if tagged)
    if mode.prot_deg_flag
        degradationRate = ...
            [tube.UserData.ReactionConfig.Protein_ClpXP_Forward tube.UserData.ReactionConfig.Protein_ClpXP_Reverse...
            tube.UserData.ReactionConfig.Protein_ClpXP_complex_deg];
        txtl_protein_degradation(mode, tube, protein,degradationRate);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_add_dna:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end


end % end of function


function [binVariable,indexes] = checkForStringInACellList(cellList,matchStr)
FlagVector = cellfun(@(x) strcmp(x,matchStr),cellList,'UniformOutput',false);
indexes = find(cell2mat(FlagVector) > 0);
if sum(cell2mat(FlagVector)) >= 1
    binVariable = true;
else
    binVariable = false;
end
end

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
