function varargout = txtl_tx_cascade(mode, tube, dna, rna, RNAP, RNAPbound, prom_spec, rbs_spec, gene_spec, varargin)
% txtl_tx_cascade.m - RNA circuit set-up file
% VS 7-25-2013
%
% Written by Richard Murray, Sep 2012
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

% cant do this because it fails in reaction mode.
% prom_spec = tube.userdata.DNAinfo{end}{1};
% rbs_spec = tube.userdata.DNAinfo{end}{2};
% gene_spec = tube.userdata.DNAinfo{end}{3};

% get the various strings specified in the promoter, utr, gene
[promData, promStr] = txtl_parsespec(prom_spec);
[utrData, utrStr] = txtl_parsespec(rbs_spec);
[geneData, geneStr] = txtl_parsespec(gene_spec);
paramObj_RNA = txtl_RNA_config('txtl_RNA');
att = utrData{1,1};
if ~strcmp(att, {'att1','att2'}) %can be extended to any number of att's
    error('Something went wrong. We think there is Attenuator RNA present');
end
% [BooleanAtt1Present,indexesOfAtt1Present] = checkForStringInACellList(utrData(1,:),'att1');
% [BooleanAnti1Present,indexesOfAnti1Present] = checkForStringInACellList(utrData(1,:),'anti1');
% SIMmoduleException = false;
% doubleAntisense = false;
% if length(indexesOfAtt1Present)==2
%     if indexesOfAtt1Present == [1 2]
%         SIMmoduleException = true;
%     end
% end
% if length(indexesOfAnti1Present)==2
%     if indexesOfAnti1Present == [2 3]
%         doubleAntisense = true;
%     end
% end

%% Setup Species
if strcmp(mode.add_dna_driver, 'Setup Species')
    RNAPbound_term = ['term_' RNAPbound];
    if nargin == 9
        if mode.double_antisense
            coreSpecies = {'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term, ['RNA ' att],['RNA ' att '-anti1'],['RNA ' att '-anti1-anti1'], RNAP};
        elseif mode.sim_module_exception
            coreSpecies = {'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term, 'RNA att1_SIM','RNA att1-att1', RNAP};
        else
            coreSpecies = {'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term, ['RNA ' att], RNAP}; 
            %! TODO 11/13/13 what sets up RNA anti? is it something in
            %txtl_add_dna? so normal transcription takes care of it? well
            %now we want this one to handle that case. or do we?
        end
    elseif nargin == 10
        extraSpecies=varargin{1};
        coreSpecies = [{'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term, ['RNA ' att], RNAP},extraSpecies]; 
    else
        error('the number of argument should be at 9 ot 10, not %d',nargin);
    end
    varargout{1} = coreSpecies;
    %% Setup Reactions
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    % there are a bunch of exceptions here. We need to generalize in a future
    % release.
    
    %if there are two sources of RNA att1, like DNA prom--att1-rbs--rfp and
    %DNA prom--att1-att1-rbs--gfp, then the two RNA att1 from these two
    %will be indistinguishable. they need to be distinguishable. 
    RNAPbound_term = ['term_' RNAPbound];
    if nargin ==9
        if mode.double_antisense % 
            transcriptionEq1 = ...
            ['[CUTP:AGTP:' RNAPbound '] -> [' RNAPbound ':RNA ' att ']'];
            transcriptionEq2 = ...
                ['[CUTP:AGTP:' RNAPbound ':RNA ' att '] -> '  RNAPbound_term ' + RNA ' att '-anti1 + RNA ' att '-anti1-anti1'];
            
            ktxExpression1 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.att'); % rna.UserData.att defined in txtl_add_dna
            ktxExpression2 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.remaining'); % rna.UserData.remaining defined in txtl_add_dna
            ktx1 = eval(ktxExpression1);
            ktx2 = eval(ktxExpression2);
            
        elseif mode.sim_module_exception % SIM module exception
            transcriptionEq1 = ...
            ['[CUTP:AGTP:' RNAPbound '] -> [' RNAPbound ':RNA att1_SIM]'];
            transcriptionEq2 = ...
                ['[CUTP:AGTP:' RNAPbound ':RNA att1_SIM] -> [' RNAPbound ':RNA att1-att1]'];
            transcriptionEq3 = ...
                ['[CUTP:AGTP:' RNAPbound ':RNA att1-att1] -> '  RNAPbound_term ' + ' rna.Name];
            
            ktxExpression1 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.attFirst');
            ktxExpression2 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.attSecond');
            ktxExpression3 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.remaining');
            ktx1 = eval(ktxExpression1);
            ktx2 = eval(ktxExpression2);
            ktx3 = eval(ktxExpression3);
            
        else % NOMINAL
            transcriptionEq1 = ...
            ['[CUTP:AGTP:' RNAPbound '] -> [' RNAPbound ':RNA ' att ']'];
            transcriptionEq2 = ...
                ['[CUTP:AGTP:' RNAPbound ':RNA ' att '] -> '  RNAPbound_term ' + ' rna.Name];
            ktxExpression1 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.att'); % rna.UserData.att defined in txtl_add_dna
            ktxExpression2 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.remaining'); % rna.UserData.remaining defined in txtl_add_dna
            ktx1 = eval(ktxExpression1);
            ktx2 = eval(ktxExpression2);
            txtl_addreaction(tube,['[' RNAPbound_term '] -> ' RNAP  ' + ' dna.Name],...
            'MassAction',{'TXTL_RNAPBOUND_TERMINATION_RATE', tube.UserData.ReactionConfig.RNAPbound_termination_rate});
        end
    elseif nargin==10
        extraSpecies = varargin{1};
            % processing the extraSpecies, no SIM modle or double antisense
            % for now. we will rewrite all this in a more general form
            % soon. 
            extraStr = extraSpecies{1};
            for k=2:size(extraSpecies,2)
                extraStr = [extraStr '+' extraSpecies{k}];
            end
        transcriptionEq1 = ...
            ['[CUTP:AGTP:' RNAPbound '] -> [' RNAPbound ':RNA ' att ']'];
            transcriptionEq2 = ...
                ['[CUTP:AGTP:' RNAPbound ':RNA ' att '] -> '  RNAPbound_term ' + ' rna.Name];
            ktxExpression1 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.att'); % rna.UserData.att defined in txtl_add_dna
            ktxExpression2 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
                'RNA_Length','rna.UserData.remaining'); % rna.UserData.remaining defined in txtl_add_dna
            ktx1 = eval(ktxExpression1);
            ktx2 = eval(ktxExpression2);
            txtl_addreaction(tube,['[' RNAPbound_term '] -> ' RNAP  ' + ' dna.Name ' + ' extraStr],...
            'MassAction',{'TXTL_RNAPBOUND_TERMINATION_RATE', tube.UserData.ReactionConfig.RNAPbound_termination_rate}); % change this rate to correct variable
    else
        error('the number of arguments should be at 9 or 10, not %d',nargin);
    end
    
    
    % binding NTP
    NTPparameters = {'TXTL_NTP_RNAP_F', tube.UserData.ReactionConfig.NTP_Forward;
        'TXTL_NTP_RNAP_R', tube.UserData.ReactionConfig.NTP_Reverse};
    NTPparameters_fast = {'TXTL_NTP_RNAP_F', 1000*tube.UserData.ReactionConfig.NTP_Forward;
        'TXTL_NTP_RNAP_R', 1000*tube.UserData.ReactionConfig.NTP_Reverse};
    
    if mode.sim_module_exception  % SIM module exception
        txtl_addreaction(tube,['[' RNAPbound '] + AGTP <-> [AGTP:' RNAPbound ']'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[' RNAPbound '] + CUTP <-> [CUTP:' RNAPbound ']'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[AGTP:' RNAPbound '] + CUTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters);
        txtl_addreaction(tube,['[CUTP:' RNAPbound '] + AGTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters);  
    
        txtl_addreaction(tube,['[' RNAPbound ':RNA att1_SIM] + AGTP <-> [AGTP:' RNAPbound ':RNA att1_SIM]'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[' RNAPbound ':RNA att1_SIM] + CUTP <-> [CUTP:' RNAPbound ':RNA att1_SIM]'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[AGTP:' RNAPbound ':RNA att1_SIM] + CUTP <-> [CUTP:AGTP:' RNAPbound ':RNA att1_SIM]'],...
            'MassAction',NTPparameters);
        txtl_addreaction(tube,['[CUTP:' RNAPbound ':RNA att1_SIM] + AGTP <-> [CUTP:AGTP:' RNAPbound ':RNA att1_SIM]'],...
            'MassAction',NTPparameters);
        
        txtl_addreaction(tube,['[' RNAPbound ':RNA att1-att1] + AGTP <-> [AGTP:' RNAPbound ':RNA att1-att1]'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[' RNAPbound ':RNA att1-att1] + CUTP <-> [CUTP:' RNAPbound ':RNA att1-att1]'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[AGTP:' RNAPbound ':RNA att1-att1] + CUTP <-> [CUTP:AGTP:' RNAPbound ':RNA att1-att1]'],...
            'MassAction',NTPparameters);
        txtl_addreaction(tube,['[CUTP:' RNAPbound ':RNA att1-att1] + AGTP <-> [CUTP:AGTP:' RNAPbound ':RNA att1-att1]'],...
            'MassAction',NTPparameters);
        
    else
        txtl_addreaction(tube,['[' RNAPbound '] + AGTP <-> [AGTP:' RNAPbound ']'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[' RNAPbound '] + CUTP <-> [CUTP:' RNAPbound ']'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[AGTP:' RNAPbound '] + CUTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters);
        txtl_addreaction(tube,['[CUTP:' RNAPbound '] + AGTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters); 
    
        txtl_addreaction(tube,['[' RNAPbound ':RNA ' att '] + AGTP <-> [AGTP:' RNAPbound ':RNA ' att ']'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[' RNAPbound ':RNA ' att '] + CUTP <-> [CUTP:' RNAPbound ':RNA ' att ']'],...
            'MassAction',NTPparameters_fast);
        txtl_addreaction(tube,['[AGTP:' RNAPbound ':RNA ' att '] + CUTP <-> [CUTP:AGTP:' RNAPbound ':RNA ' att ']'],...
        'MassAction',NTPparameters);
        txtl_addreaction(tube,['[CUTP:' RNAPbound ':RNA ' att '] + AGTP <-> [CUTP:AGTP:' RNAPbound ':RNA ' att ']'],...
        'MassAction',NTPparameters); 
    end
    
    
    % dummy reaction for NTP consumption
    if mode.sim_module_exception  % SIM module exception
        ntpcnt1 = ceil(rna.UserData.attFirst/2);   	% get number of NTP blocks
        ntpcnt2 = ceil(rna.UserData.attSecond/2);   	% get number of NTP blocks
        ntpcnt3 = ceil(rna.UserData.remaining/2);	% get number of NTP blocks
        NTPConsumptionRate1 = {'TXTL_NTP_consumption1',(ntpcnt1-1)*ktx1};
        NTPConsumptionRate2 = {'TXTL_NTP_consumption2',(ntpcnt2-1)*ktx2};
        NTPConsumptionRate3 = {'TXTL_NTP_consumption2',(ntpcnt3-1)*ktx3};
        txtl_addreaction(tube,['[CUTP:AGTP:' RNAPbound '] -> ' RNAPbound],...
            'MassAction',NTPConsumptionRate1);
        txtl_addreaction(tube,['[CUTP:AGTP:' RNAPbound ':RNA att1_SIM] -> [' RNAPbound ':RNA att1_SIM]'],...
            'MassAction',NTPConsumptionRate2);
        txtl_addreaction(tube,['[CUTP:AGTP:' RNAPbound ':RNA att1-att1] -> [' RNAPbound ':RNA att1_att1]'],...
            'MassAction',NTPConsumptionRate3);
    else
        ntpcnt1 = ceil(rna.UserData.att/2);   	% get number of NTP blocks
        ntpcnt2 = ceil(rna.UserData.remaining/2);	% get number of NTP blocks
        NTPConsumptionRate1 = {'TXTL_NTP_consumption1',(ntpcnt1-1)*ktx1};
        NTPConsumptionRate2 = {'TXTL_NTP_consumption2',(ntpcnt2-1)*ktx2};
        txtl_addreaction(tube,['[CUTP:AGTP:' RNAPbound '] -> ' RNAPbound],...
            'MassAction',NTPConsumptionRate1);
        txtl_addreaction(tube,['[CUTP:AGTP:' RNAPbound ':RNA ' att '] -> [' RNAPbound ':RNA ' att ']' ],...
            'MassAction',NTPConsumptionRate2);
    end
    
    %att anti reactions
    attnumber = att(end);
    [~,listOfSpecies] = getstoichmatrix(tube);
    matchStr = regexp(listOfSpecies,['(^RNA .*anti' attnumber '.*)'],'tokens','once'); % this end-anchoring of anti must be removed to allow for anti to be in the middle
    tempList = vertcat(matchStr{:});
     matchStr2 = regexp(tempList,':');
    antiIdx = cellfun(@isempty, matchStr2);
    listOfantisense = tempList(antiIdx);
    attanti = ['att' attnumber '_anti' attnumber];
    
    if mode.sim_module_exception
        if ~isempty(listOfantisense)
            for k = 1:size(listOfantisense,1)
                if strcmp(listOfantisense{k}, 'RNA att2-anti1-anti1') %|| strcmp(listOfantisense{k}, 'RNA att2-anti1')
                    complex1_rate = {'TXTL_RNA_C1_F', paramObj_RNA.att1_anti11_repression_F_rate;'TXTL_RNA_C1_R', paramObj_RNA.att1_anti11_repression_R_rate};
                    att_anti_termination_rate = {'TXTL_RNA_ATTANTI_TERM', paramObj_RNA.att1_anti11_termination_rate};
                else
                    eval(['complex1_rate = {''TXTL_RNA_C1_F'', paramObj_RNA.' attanti '_repression_F_rate;''TXTL_RNA_C1_R'', paramObj_RNA.' attanti '_repression_R_rate};'])
                    eval(['att_anti_termination_rate = {''TXTL_RNA_ATTANTI_TERM'', paramObj_RNA.' attanti '_termination_rate};'])
                end                % first repression
                txtl_addreaction(tube,...
                    [RNAPbound ':RNA att1_SIM + ' listOfantisense{k} ' <-> [' RNAPbound ':RNA att1_SIM:' listOfantisense{k} ']'],...
                    'MassAction',complex1_rate);
                txtl_addreaction(tube,['[' RNAPbound ':RNA att1_SIM:' listOfantisense{k} '] -> ' dna.Name ' + [' listOfantisense{k} ':RNA att1_SIM] + ' RNAP],...
                    'MassAction',att_anti_termination_rate);
                txtl_addreaction(tube,...
                    [RNAPbound ':RNA att1-att1 + ' listOfantisense{k} ' <-> [' RNAPbound ':RNA att1-att1:' listOfantisense{k} ']'],...
                    'MassAction',complex1_rate);
                txtl_addreaction(tube,['[' RNAPbound ':RNA att1-att1:' listOfantisense{k} '] -> ' dna.Name ' + [' listOfantisense{k} ':RNA att1-att1] + ' RNAP],...
                    'MassAction',att_anti_termination_rate);
                
                for j = 1:size(listOfantisense,1)
                    if strcmp(listOfantisense{j}, 'RNA att2-anti1-anti1') %|| strcmp(listOfantisense{j}, 'RNA att2-anti1')
                        complex1_rate = {'TXTL_RNA_C1_F', paramObj_RNA.att1_anti11_repression_F_rate;'TXTL_RNA_C1_R', paramObj_RNA.att1_anti11_repression_R_rate};
                        att_anti_termination_rate = {'TXTL_RNA_ATTANTI_TERM', paramObj_RNA.att1_anti11_termination_rate};
                    else
                        eval(['complex1_rate = {''TXTL_RNA_C1_F'', paramObj_RNA.' attanti '_repression_F_rate;''TXTL_RNA_C1_R'', paramObj_RNA.' attanti '_repression_R_rate};'])
                        eval(['att_anti_termination_rate = {''TXTL_RNA_ATTANTI_TERM'', paramObj_RNA.' attanti '_termination_rate};'])
                    end
                    % second repression
                    txtl_addreaction(tube,...
                        [RNAPbound ':RNA att1-att1:' listOfantisense{k} ' + ' listOfantisense{j} ' <-> [' RNAPbound ':RNA att1-att1:' listOfantisense{k} ':' listOfantisense{j} ']'],...
                        'MassAction',complex1_rate);
                    txtl_addreaction(tube,['[' RNAPbound ':RNA att1-att1:' listOfantisense{k} ':' listOfantisense{j} '] -> ' dna.Name ' + [' listOfantisense{k} ':' listOfantisense{j} ':RNA att1-att1] + ' RNAP],...
                        'MassAction',att_anti_termination_rate);
                end
                
            end
        end
        
    else
        
        if ~isempty(listOfantisense)
            for k = 1:size(listOfantisense,1)
                if strcmp(listOfantisense{k}, 'RNA att2-anti1-anti1') %|| strcmp(listOfantisense{k}, 'RNA att2-anti1')
                    complex1_rate = {'TXTL_RNA_C1_F', paramObj_RNA.att1_anti11_repression_F_rate;'TXTL_RNA_C1_R', paramObj_RNA.att1_anti11_repression_R_rate};
                    att_anti_termination_rate = {'TXTL_RNA_ATTANTI_TERM', paramObj_RNA.att1_anti11_termination_rate};
                else
                    eval(['complex1_rate = {''TXTL_RNA_C1_F'', paramObj_RNA.' attanti '_repression_F_rate;''TXTL_RNA_C1_R'', paramObj_RNA.' attanti '_repression_R_rate};'])
                    eval(['att_anti_termination_rate = {''TXTL_RNA_ATTANTI_TERM'', paramObj_RNA.' attanti '_termination_rate};'])
                end
                txtl_addreaction(tube,...
                    [RNAPbound ':RNA ' att ' + ' listOfantisense{k} ' <-> [' RNAPbound ':RNA ' att ':' listOfantisense{k} ']'],...
                    'MassAction',complex1_rate);
                txtl_addreaction(tube,['[' RNAPbound ':RNA ' att ':' listOfantisense{k} '] -> ' dna.Name ' + [' listOfantisense{k} ':RNA ' att '] + ' RNAP],...
                    'MassAction',att_anti_termination_rate);
            end
        end
    end
    
    
    
%     %att anti cross talk reactions
%     if attnumber == '1'
%         attother = '2';
%     elseif attnumber == '2'
%         attother = '1';
%     else
%         error('txtltoolbox:txtl_tx_cascade:undefinedatt', ...
%             'The possible attenuators are att1 and att2');
%     end
%     
%     [~,listOfSpecies] = getstoichmatrix(tube);
%     matchStr = regexp(listOfSpecies,['(^RNA .*anti' attother '$)'],'tokens','once');
%     listOfantisense = vertcat(matchStr{:});
%     
%     
%     if mode.sim_module_exception
%         if ~isempty(listOfantisense)
%             for k = 1:size(listOfantisense,1)
%                 if strcmp(listOfantisense{k}, 'RNA att2-anti1-anti1') %|| strcmp(listOfantisense{k}, 'RNA att2-anti1')
%                     mRNA_Xtalk_rate = {'TXTL_RNA_XTALK_F', paramObj_RNA.att2_anti1_anti1_Xtalk_F_rate;'TXTL_RNA_XTALK_R',paramObj_RNA.att2_anti1_anti1_Xtalk_R_rate};
%                     att_anti_Xtalk_termination_rate = {'TXTL_RNA_ATTANTI_XTALK_TERM', paramObj_RNA.att2_anti11_Xtalk_termination_rate};
%                 else
%                     eval(['mRNA_Xtalk_rate = {''TXTL_RNA_XTALK_F'', paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_F_rate;''TXTL_RNA_XTALK_R'',paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_R_rate};'])
%                     eval(['att_anti_Xtalk_termination_rate = {''TXTL_RNA_ATTANTI_XTALK_TERM'', paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_termination_rate};'])
%                 end
%                 txtl_addreaction(tube,...
%                     [RNAPbound ':RNA att1_SIM + ' listOfantisense{k} ' <-> [' RNAPbound ':RNA att1_SIM:' listOfantisense{k} ']'],...
%                     'MassAction',mRNA_Xtalk_rate);
%                 txtl_addreaction(tube,['[' RNAPbound ':RNA att1_SIM:' listOfantisense{k} '] -> ' dna.Name ' + [' listOfantisense{k} ':RNA att1_SIM] + ' RNAP],...
%                     'MassAction',att_anti_Xtalk_termination_rate);
%                 txtl_addreaction(tube,...
%                     [RNAPbound ':RNA att1-att1 + ' listOfantisense{k} ' <-> [' RNAPbound ':RNA att1-att1:' listOfantisense{k} ']'],...
%                     'MassAction',mRNA_Xtalk_rate);
%                 txtl_addreaction(tube,['[' RNAPbound ':RNA att1-att1:' listOfantisense{k} '] -> ' dna.Name ' + [' listOfantisense{k} ':RNA att1-att1] + ' RNAP],...
%                     'MassAction',att_anti_Xtalk_termination_rate);
%                 
%                 for j = 1:size(listOfantisense,1)
%                     % SIM module exception
%                     if strcmp(listOfantisense{k}, 'RNA att2-anti1-anti1') %|| strcmp(listOfantisense{k}, 'RNA att2-anti1')
%                         mRNA_Xtalk_rate = {'TXTL_RNA_XTALK_F', paramObj_RNA.att2_anti1_anti1_Xtalk_F_rate;'TXTL_RNA_XTALK_R',paramObj_RNA.att2_anti1_anti1_Xtalk_R_rate};
%                         att_anti_Xtalk_termination_rate = {'TXTL_RNA_ATTANTI_XTALK_TERM', paramObj_RNA.att2_anti11_Xtalk_termination_rate};
%                     else
%                         eval(['mRNA_Xtalk_rate = {''TXTL_RNA_XTALK_F'', paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_F_rate;''TXTL_RNA_XTALK_R'',paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_R_rate};'])
%                         eval(['att_anti_Xtalk_termination_rate = {''TXTL_RNA_ATTANTI_XTALK_TERM'', paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_termination_rate};'])
%                     end
%                     txtl_addreaction(tube,...
%                         [RNAPbound ':RNA att1-att1:' listOfantisense{k} ' + ' listOfantisense{j} ' <-> [' RNAPbound ':RNA att1-att1:' listOfantisense{k} ':' listOfantisense{j} ']'],...
%                         'MassAction',mRNA_Xtalk_rate);
%                     txtl_addreaction(tube,['[' RNAPbound ':RNA att1-att1:' listOfantisense{k} ':' listOfantisense{j} '] -> ' dna.Name ' + [' listOfantisense{k} ':' listOfantisense{j} ':RNA att1-att1] + ' RNAP],...
%                         'MassAction',att_anti_Xtalk_termination_rate);
%                 end
%                 
%             end
%         end
%     else
%         if ~isempty(listOfantisense)
%             for k = 1:size(listOfantisense,1)
%                 if strcmp(listOfantisense{k}, 'RNA att2-anti1-anti1') %|| strcmp(listOfantisense{k}, 'RNA att2-anti1')
%                     mRNA_Xtalk_rate = {'TXTL_RNA_XTALK_F', paramObj_RNA.att2_anti1_anti1_Xtalk_F_rate;'TXTL_RNA_XTALK_R',paramObj_RNA.att2_anti1_anti1_Xtalk_R_rate};
%                     att_anti_Xtalk_termination_rate = {'TXTL_RNA_ATTANTI_XTALK_TERM', paramObj_RNA.att2_anti11_Xtalk_termination_rate};
%                 else
%                     eval(['mRNA_Xtalk_rate = {''TXTL_RNA_XTALK_F'', paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_F_rate;''TXTL_RNA_XTALK_R'',paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_R_rate};'])
%                     eval(['att_anti_Xtalk_termination_rate = {''TXTL_RNA_ATTANTI_XTALK_TERM'', paramObj_RNA.att' attnumber '_anti' attother '_Xtalk_termination_rate};'])
%                 end
%                 txtl_addreaction(tube,...
%                     [RNAPbound ':RNA ' att ' + ' listOfantisense{k} ' <-> [' RNAPbound ':RNA ' att ':' listOfantisense{k} ']'],...
%                     'MassAction',mRNA_Xtalk_rate);
%                 txtl_addreaction(tube,['[' RNAPbound ':RNA ' att ':' listOfantisense{k} '] -> ' dna.Name ' + [' listOfantisense{k} ':RNA ' att '] + ' RNAP],...
%                     'MassAction',att_anti_Xtalk_termination_rate);
%                 
%             end
%         end
%     end
    % Auto-Termination att1_auto_termination_rate
    eval(['auto_termination_rate = {''TXTL_RNA_ATT_AUTOTERM'',paramObj_RNA.att' attnumber '_auto_termination_rate};'])

    if mode.sim_module_exception  % SIM module exception
        % transcription
        txtl_addreaction(tube,transcriptionEq1,'MassAction',{'TXTL_transcription_rate1',ktx1});
        txtl_addreaction(tube,transcriptionEq2,'MassAction',{'TXTL_transcription_rate2',ktx2});
        txtl_addreaction(tube,transcriptionEq3,'MassAction',{'TXTL_transcription_rate2',ktx3});
        % Auto-Termination att1_auto_termination_rate
        txtl_addreaction(tube,['[' RNAPbound ':RNA att1_SIM] -> '  RNAPbound ' + [RNA att1_SIM]'],...
            'MassAction',auto_termination_rate);
        txtl_addreaction(tube,['[' RNAPbound ':RNA att1-att1] -> '  RNAPbound ' + [RNA att1-att1]'],...
            'MassAction',auto_termination_rate);
        
    else
        % transcription
        txtl_addreaction(tube,transcriptionEq1,'MassAction',{'TXTL_transcription_rate1',ktx1});
        txtl_addreaction(tube,transcriptionEq2,'MassAction',{'TXTL_transcription_rate2',ktx2});
        txtl_addreaction(tube,['[' RNAPbound ':RNA ' att '] -> [RNA ' att '] + ' RNAPbound],...
            'MassAction',auto_termination_rate);
    end
    
else
    error('txtltoolbox:txtl_tx_cascade:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end
end

%
% function [binVariable,indexes] = checkForStringInACellList(cellList,matchStr)
% FlagVector = cellfun(@(x) strcmp(x,matchStr),cellList,'UniformOutput',false);
% indexes = find(cell2mat(FlagVector) > 0);
% if sum(cell2mat(FlagVector)) >= 1
%     binVariable = true;
% else
%     binVariable = false;
% end
% end
