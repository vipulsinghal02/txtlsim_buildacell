%%
% This function checks that every protein has its the corresponding species/reactions
% in the species/reaction list. For example, if one manually adds a speices to the
% TXTL reaction all the corresponding species/reactions will be set up in the is
% function
%
function txtl_setup_new_protein_added(tube, add_dna_mode)

% Step 1: Make a list of all the proteins added by DNA strings:
speciesNames = get(tube.species, 'name');
matchStr = regexp(speciesNames,'(^DNA (\w*(--|-)?)*)','tokens','once');
listOfDNA = vertcat(matchStr{:});
listOfDNA = unique(listOfDNA);
DNAparts = regexp(listOfDNA, '--','split');
% If a DNA constract has gene on it, it will have 3 parts, let's select
% those
DNAwithGene = cellfun(@(x) size(x,2)>2, DNAparts);
DNAparts = DNAparts(DNAwithGene > 0);
DNAparts = vertcat(DNAparts{:});
proteinsAlreadySetUp = DNAparts(:,3);

% Step2: Make a list of protein in the model (based on the species list):
matchStr = regexp(speciesNames,'(^protein.*)','tokens','once'); 
listOfprotein = vertcat(matchStr{:});

% Step3: Compare the list of proteins and the list of genes from the DNA
% string
matchStr = regexp(listOfprotein,'^protein (\w*(-lva)?(-terminator)?)','tokens','once');
matchStr = vertcat(matchStr{:});
proteinList = regexprep(matchStr,'tetramer|dimer','');
proteinList = unique(proteinList);

%implementation of union has changed in 2013a [~, ia, ~] = union(proteinList,proteinsAlreadySetUp);
ia = [];
for k=1:size(proteinList,1)
    if findStringInAList(proteinsAlreadySetUp,proteinList(k)) == 0
        ia(end+1) = k;
    end
end


% Proteins from the first list have not been set up (indexed by ia)
for k = 1:size(ia,1)
    
    fileName = strrep(proteinList{ia(k)},'-lva','');
    if exist(['txtl_protein_'  fileName], 'file') == 2
        % Run the protein specific setup
        protData = txtl_parsespec(proteinList{ia(k)});
        proteinIdx = findspecies(tube, ['protein ' proteinList{ia(k)}]);
        % handling the case when the protein is not exist in the species list, e.g. a dimer was
        % added then the monomer is not exist yet.
        if proteinIdx == 0
            proteinObj = txtl_addspecies(tube,['protein ' proteinList{ia(k)}],0);
        else
            proteinObj = tube.Species(proteinIdx);
        end
        
        eval(['txtl_protein_' fileName '(add_dna_mode, tube, proteinObj, protData);']);
        
        % double check for degradation tags, if necessary set up degradation
        %
        if ~isempty(strfind(proteinList{ia(k)},'-lva'))
            degradationRate = ...
                [tube.UserData.ReactionConfig.Protein_ClpXP_Forward tube.UserData.ReactionConfig.Protein_ClpXP_Reverse...
                tube.UserData.ReactionConfig.Protein_ClpXP_complex_deg];
            
            txtl_protein_degradation(add_dna_mode, tube, proteinObj,degradationRate);
        end
        
    end
    
end

end