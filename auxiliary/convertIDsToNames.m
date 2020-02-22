function convertIDsToNames(translatedSBML,writtenODEfile,varargin)

% 
if nargin > 2 
    outputFileName = varargin{1};
else
    outputFileName = '%s_odefile.m';
end
%% build associative map
idVsName = {};

% compartments
idVsName = vertcat(idVsName,getIdNamePairs(translatedSBML.compartment));

% species
idVsName = vertcat(idVsName,getIdNamePairs(translatedSBML.species));

% parameters
NumOfReactions = size(translatedSBML.reaction,2);
for k=1:NumOfReactions
    idVsName = vertcat(idVsName,getIdNamePairs(translatedSBML.reaction(k).kineticLaw.parameter));
end
% reaction ids 
replaceIDs = cellfun(@num2str,mat2cell([1:19]',[ones(19,1)],1),'UniformOutput',false);
tmp = getIdNamePairs(translatedSBML.reaction);
tmp{2} = replaceIDs;
idVsName = vertcat(idVsName,tmp);


% modell name
idVsName = vertcat(idVsName,getIdNamePairs(translatedSBML));
fileName = idVsName{end,2}{1};



ids = vertcat(idVsName{:,1});
name = vertcat(idVsName{:,2}) ;

%% readfile line-by-line and search for ids to replace

if exist(writtenODEfile,'file')
    fid = fopen(writtenODEfile,'r');
    fid2 = fopen(sprintf(outputFileName,fileName),'w');
    while ~feof(fid)
        s = fgetl(fid);
        %z = strrep(s,ids,name)
        
        z  = regexprep(s,ids,name);
        fprintf(fid2,'%s\n',z);
        
    end
    
    fclose(fid);
    fclose(fid2);
else
    error('no such file: %s',writtenODEfile);
end



%% aux function
    function pairs = getIdNamePairs(structWithData)
        NumOfItems = size(structWithData,2);
        tmp = cell(NumOfItems,2);
        for l = 1:NumOfItems
            tmp{l,1} = structWithData(l).id;
            tmp{l,2} = regexprep(structWithData(l).name,{' ',':','--','*','-'},{'_'});
        end
        pairs = {vertcat(tmp(:,1)) vertcat(tmp(:,2))};
    end
end



