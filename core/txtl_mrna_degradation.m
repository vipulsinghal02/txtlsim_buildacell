function txtl_mrna_degradation(mode, tube, dna, rna, rbs_spec)
% Vipul Singhal Feb 13 2014
% We model RNA degradation as an enzymatic step reaction, as described in
% Vincent Noireaux's 2011 PRL paper 
% "Coarse-Grained Dynamics of Protein Synthesis in a Cell-Free System"
% 200 nM of radioactively labelled 960bases long RNA decause exponentially
% with a half life of 12 min. This means that the RNAse does not get
% saturated at this concentration of mRNA. 


% includeRNase = false;
% BooleanAtt1Present = false;
% BooleanAtt2Present = false;
% BooleanAnti1Present = false;
% BooleanAnti2Present = false;
% for i = 1: length(tube.userdata.DNAinfo)
%     UTR = tube.userdata.DNAinfo{i}{2};
%     [utrData, utrStr] = txtl_parsespec(UTR);
%     [BooleanAtt1Present,indexesOfAtt1Present] = checkForStringInACellList(utrData(1,:),'att1');
%     [BooleanAnti1Present,indexesOfAnti1Present] = checkForStringInACellList(utrData(1,:),'anti1');
%         [BooleanAtt2Present,indexesOfAtt2Present] = checkForStringInACellList(utrData(1,:),'att2');
%     [BooleanAnti2Present,indexesOfAnti2Present] = checkForStringInACellList(utrData(1,:),'anti2');
% if BooleanAtt1Present || BooleanAtt2Present || BooleanAnti1Present || BooleanAnti2Present
%     includeRNase = true;
% end
% end


    complexF = tube.UserData.ReactionConfig.RNase_F;
    complexR = tube.UserData.ReactionConfig.RNase_R;
    degRate = tube.UserData.ReactionConfig.RNA_deg;
    if isempty(tube.UserData.ReactionConfig.RNase_F)
        complexF = degRate/10;
        complexR = degRate/40;
    end
    length_over_4 = round(rna.UserData/4);
    % Setup RNA degradation reactions, by searching for strings with rna in
    % them, and then degrading those strings. 
    listOfSpecies = get(tube.species, 'name');
    rnastr = ['(' rna.Name '(?=:))|(' rna.Name ')'] ; 
    % only need the second one actually, dont need the lookahead for the colon
    
    rnalist = regexp(listOfSpecies, rnastr, 'match');
    rna_species_ix = ~cellfun(@isempty,rnalist);
    RNAcomplexes = listOfSpecies(rna_species_ix);
    [splset] = regexp(RNAcomplexes,  '\:', 'split');
    
    for i = 1:length(splset)
        fullsplit = splset{i}';
        rnasplit = regexp(fullsplit, rnastr, 'match');
        nonRNA_ix = cellfun(@isempty,rnasplit);
        RNA_ix = ~cellfun(@isempty,rnasplit);
        txtl_addreaction(tube,[RNAcomplexes{i} ' + RNase <-> ' RNAcomplexes{i} ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        nonRNAlist = fullsplit(nonRNA_ix);
        
        productspecies = [];
        for j = 1:length(nonRNAlist)
            % if in the complex we have term_Ribo, then we need to
            % return Ribo, not create nonsensical species term_Ribo. 
            if strcmp(nonRNAlist{j}, 'term_Ribo')
                nonRNAlist{j} = 'Ribo';
            end
            productspecies = [productspecies ' + ' nonRNAlist{j}];
        end
        
        if isfield(tube.UserData, 'energymode') && strcmp(tube.UserData.energymode, 'regeneration')
           
                txtl_addreaction(tube,[RNAcomplexes{i} ':RNase -> RNase + '...
                    num2str(length_over_4) ' AGMP + ' num2str(length_over_4) ' CUMP '...
                    productspecies],...
                    'MassAction',{'TXTL_RNAdeg_kc',degRate});
            
        else
            txtl_addreaction(tube,[RNAcomplexes{i} ':RNase -> RNase ' productspecies],...
                'MassAction',{'TXTL_RNAdeg_kc',degRate});
        end
    end
end

function [binVariable,indexes] = checkForStringInACellList(cellList,matchStr)
FlagVector = cellfun(@(x) strcmp(x,matchStr),cellList,'UniformOutput',false);
indexes = find(cell2mat(FlagVector) > 0);
if sum(cell2mat(FlagVector)) >= 1
    binVariable = true;
else
    binVariable = false;
end
end