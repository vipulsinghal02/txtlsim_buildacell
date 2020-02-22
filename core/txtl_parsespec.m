%% Utility function for parsing out a specification string
function [parsedData, combinedStr] = txtl_parsespec(spec)

indivRegions = regexp(spec, '-','split'); %cell array of individual xyz(123) strings
namesAndLengths = regexp(indivRegions, '\w+','match'); %cell array of cells containing the names and lengths of dna regions
names = cell(1,length(namesAndLengths));
lengths = cell(1,length(namesAndLengths));
combinedStr = '';

%error checking followed by returning parsed strings
for k = 1:length(namesAndLengths)
    if isempty(namesAndLengths{k})
        error('txtl_add_dna:wrongStringFormat',...
            ['the string %s should be: name(length)-name2(length2)-' ...
            '...-nameN(lengthN), where the lengths are optional. eg: thio-junk(500)-ptet(50)'...
            'the name must start with an alphabet'], spec)
    else
        A = isstrprop(namesAndLengths{k}{1},'alpha');
        if ~A(1) % this happens when the name of the dna fragment does not start with an alphabet
            error('txtl_add_dna:wrongSpeciesName',...
                ['species named %s should start with an alphabet. Format is' ...
                ' name(length). Where the lengths are optional. eg: thio or junk(500)'],indivRegions{k})
        end
    end
    % return the parsed name and optional length
    names{k} = namesAndLengths{k}{1};
    if length(namesAndLengths{k}) == 1
        lengths{k} = [];
    else if length(namesAndLengths{k})>2
            error('txtl_add_dna:tooManyElements',...
                ['the string %s is not of the format name(length). '...
                'It has unwanted elements after '')'''],...
                indivRegions{k});
        else if length(namesAndLengths{k})==2
                lengths{k} = str2double(namesAndLengths{k}{2});
            end
        end
    end
    if k==1
        combinedStr = names{k};
    else
        combinedStr = [combinedStr '-' names{k}];
    end
    parsedData = [names;lengths];
end
% !TODO add error checking for numerical values for the lengths.

end