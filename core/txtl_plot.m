function processedData = txtl_plot(varargin)
% initial version for txtl_plot, the RNAs and proteins are automatically
% exploited from the provided DNA sequence.
% t_ode: nx1 time vector, no time scaling is applied inside!
% x_ode: nxm species vector
% modelObj: simBiology object of the current model
% dataGroups: special data structure for the plots 
% (regular expressions can be added)
%
% example of the data structure
%  * plotName is the name of desired plot (available plots are listed below)
%
%   Currently 3 types of plots are supported:
%   * DNA and mRNA plot (case sensitive!)
%   * Gene Expression plot
%   * Resource usage
%
%  * SpeciesToPlot: name-list of the simulation data, which will be plotted.
%    (RNA and protein names are exploited and plotted automatically from DNA sequences)
%  +-----------------------------------------
%  |Special strings:
%  |
%  | * it handles keywords (e.g. ALL_DNA/ALL_PROTEIN: plots all available dns/protein species)
%  | * it alse handles matlab compatible regular expressions
%  |   (e.g. plotting all of the protein in the system: dataGroups{2,2} = {'#(protein \w*)'};)
%  | * txtl_plot also can calculate the total concentration of selected
%  | proteins and its variants with sprint "[protein name]_tot", where protein
%  | is a valid Species name in the modelObj. (e.g. dataGroups{2,2} = {'[protein lacI]_tot'} )
%  |
%  +-----------------------------------------
%
%  * colorCodes is designated for user defined line style and coloring
%

%

%%

% get default data structrure for plot
defaultdataGroups = txtl_getDefaultPlotDataStruct();

operationMode = 'standalone';

switch nargin
    
    case 2
        simData = varargin{1};
        modelObj = varargin{2};
        t_ode = simData.Time;
        x_ode = simData.Data;
        dataGroups = defaultdataGroups;
    case 3
        if isa(varargin{1},'SimData')
            simData = varargin{1};
            t_ode = simData.Time;
            x_ode = simData.Data;
            modelObj = varargin{2};
            dataGroups = varargin{3};
        else
            t_ode = varargin{1};
            x_ode = varargin{2};
            modelObj = varargin{3};
            dataGroups = defaultdataGroups;
        end
    case 4
        t_ode = varargin{1};
        x_ode = varargin{2};
        modelObj = varargin{3};
        dataGroups = varargin{4};
    case 5
        t_ode = varargin{1};
        x_ode = varargin{2};
        modelObj = varargin{3};
        dataGroups = varargin{4};
        axesHandles = varargin{5};
        operationMode = 'GUI';
    otherwise
        error('');
end

numOfGroups = size(dataGroups,2);
listOfProteins = {};
listOfRNAs = {};
listOfDNAs = {};
if size(modelObj.Compartments,1) > 1
    currModelObj = modelObj.Compartments(1);
else 
    currModelObj = modelObj;
end

listOfSpecies = get(currModelObj.species, 'name');


if strcmp(operationMode,'standalone')
    figure('Name',currModelObj.Name); clf();
end

% building the output cell structure
processedData = cell(1,3);

% Keywords lookup table

keywords = {'ALL_DNA','#(^DNA (\w+[-=]*)*)'; 'ALL_PROTEIN','#(^protein (\w+[-=]*\*{0,3})*)'};

for k = 1:numOfGroups
    
    %%%%% DNA and mRNA plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(strcmp(dataGroups(k).plotName,'DNA and mRNA'))
        
        %! TODO further refinement of str spliting
        currentSpeciesToPlot = dataGroups(k).SpeciesToPlot;
        if ~isempty(currentSpeciesToPlot)
            % search for keywords
            [row,col] = find(cell2mat(cellfun(@(x) strcmp(currentSpeciesToPlot,x),keywords(:,1),'UniformOutput',false)) == 1);
            currentSpeciesToPlot(col) = keywords(row,2);
            regexp_ind = strmatch('#', currentSpeciesToPlot);
            if ~isempty(regexp_ind)
                autoSpecies = extractRegexpAndExecute(regexp_ind,listOfSpecies,currentSpeciesToPlot);
                % remove regexp str form dataGroups array
                currentSpeciesToPlot(regexp_ind) = [];
                listOfDNAs =  horzcat(listOfDNAs,autoSpecies);
                % delete multiple occurrences
                listOfDNAs =  unique(listOfDNAs);
            end
            listOfDNAs = horzcat(listOfDNAs,currentSpeciesToPlot);
            
            r = regexp(listOfDNAs,'--','split');
            %! TODO skip already added proteins and mRNAs to avoid duplicates
            % get each RNAs
            RNAs = cellfun(@(x) strcat('RNA',{' '},x(2),'--',x(3)),r);
            listOfRNAs = horzcat(listOfRNAs,RNAs);
            % calculate the total amount for each RNA
            totRNAs = cellfun(@(x) totalAmountOfSpecies(listOfSpecies,x_ode,x),RNAs,'UniformOutput',false);
            % merge into one cell array
            totRNAs = vertcat(totRNAs{:});
            % get each proteins
            proteins = cellfun(@(x) strcat('protein',{' '}, x(3)), r); % adding protein string for each element
            listOfProteins = horzcat(listOfProteins,proteins');
            
            
            % collect the data
            listOfDNAsRNAs = vertcat(listOfDNAs,listOfRNAs);
            dataX = getDataForSpecies(currModelObj,x_ode,listOfDNAsRNAs');
            dataDNAs = getDataForSpecies(currModelObj,x_ode,listOfDNAs');
            dataRNAs = getDataForSpecies(currModelObj,x_ode,listOfRNAs');
            
            % replace free RNA concentration with total RNA concentration
            [~, ia, ib] = intersect(listOfRNAs, totRNAs(:,1));
            dataRNAs(:,ia) = horzcat(totRNAs{ib,2});
            
        else
            warning('No DNA strings were provided!');
        end
        % plot the data
        
        % ---- Calling the txtl_plot standalone -> figure is generated -------%
        if strcmp(operationMode,'standalone')
            currentHandler = subplot(223);
            % ---- GUI mode -> graphic data is given to the appropriate handler --%
        else
            currentHandler = axesHandles.dnaRna;
        end
        
        if (~isempty(dataGroups(k).colorCodes))
            [ColorMtx,LineStyle] = getColorAndLineOrderByUserData(dataGroups(k).colorCodes);
            ax2 = axes('Position',get(currentHandler,'Position'),...
                'YAxisLocation','right',...
                'Color','none',...
                'YColor','k','XTick',[]);
            
            hold(currentHandler,'on');
            for l=1:size(dataDNAs,2)
                hl2 = line(t_ode/3600,dataRNAs(:,l),'Parent',currentHandler,'Color',ColorMtx(l,:),'LineWidth',1);
                hl1 = line(t_ode/3600,dataDNAs(:,l),'Parent',ax2,'Color',ColorMtx(l,:),'LineStyle','--');
            end
            
            hold(currentHandler,'off');
        else
            plot(currentHandler,t_ode/60,dataX);
        end
        
        lgh =legend(currentHandler,listOfRNAs, 'Location', 'NorthEast');
        if verLessThan('matlab', '8.4')
            legend(lgh,'boxoff');
        else
            lgh.Box = 'off';
        end
        
        ylabel(currentHandler,'mRNA amounts [nM]');
        ylabel(ax2,'DNA amounts [nM]');
        xlabel(currentHandler,'Time [h]');
        title(currentHandler,dataGroups(k).plotName);
        
        % add the processed data to the output structure
        processedData{2} = {['Time' listOfDNAsRNAs'],[t_ode dataX]};
        
        %%%%%%% Gene Expression plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif(strcmp(dataGroups(k).plotName,'Gene Expression'))
        
        % if total amount of selected protein is calculated
        totalAmount = {};
        % add extra user defined proteins
        currentSpeciesToPlot = dataGroups(k).SpeciesToPlot;
        if ~isempty(currentSpeciesToPlot)
            %%% search for keywords
            [row,col] = find(cell2mat(cellfun(@(x) strcmp(currentSpeciesToPlot,x),keywords(:,1),'UniformOutput',false)) == 1);
            currentSpeciesToPlot(col) = keywords(row,2);
            
           
            %%% part1 for calculating the total concentration of selected proteins
            matchStrings = regexp(currentSpeciesToPlot,'^\[(protein \w*)\]_tot','tokens','once');
            % saving the indexies of special strings
            needlessStr = find(cellfun(@(x) isempty(x),matchStrings) == 0);
            % combining the result into one cellarray;
            matchStrings = vertcat(matchStrings{:});
            % deleting special strings
            currentSpeciesToPlot(needlessStr) = [];
     
            %%% checking for special strings
            regexp_ind = strmatch('#', currentSpeciesToPlot);
            
            if ~isempty(regexp_ind)
                % if regular expression is present, than it is executed and the result goes to the listOfProteins
                auto_species =  extractRegexpAndExecute(regexp_ind,listOfSpecies,currentSpeciesToPlot);
                % remove regexp str form dataGroups array
                currentSpeciesToPlot(regexp_ind) = [];
                
                listOfProteins =  horzcat(listOfProteins,currentSpeciesToPlot,auto_species');
                % delete multiple occurrences
                listOfProteins =  unique(listOfProteins);
            end    
             
            %%% part2 calculating the total concentration of selected proteins

            if ~isempty(matchStrings)
                
                %totalAmount = cell(size(matchStrings,1),2);
                for z = 1:size(matchStrings,1)
                    tA = totalAmountOfSpecies(listOfSpecies,x_ode,matchStrings{z});
                    totalAmount = vertcat(totalAmount,tA);
                end % end for z =
                
            end
            
            %%% finally, add the manually listed species
            listOfProteins =  horzcat(listOfProteins,currentSpeciesToPlot);
            
        end
        dataX = getDataForSpecies(currModelObj,x_ode,listOfProteins);
        % adding total protein concentraion into the common data matrix and a
        % label matrix as well (This could be done before, because the
        % listOfProteins was used for aquiring Species data by name)
        if size(totalAmount,1) > 0
            for k = 1:size(totalAmount,1)
                listOfProteins(end+1) = {sprintf('[%s]_{tot}',totalAmount{k,1})};
                dataX(:,end+1) = totalAmount{k,2};
            end
        end
        
        % ---- Calling the txtl_plot standalone -> figure is generated -------%
        if strcmp(operationMode,'standalone')
            currentHandler = subplot(2,2,1:2);
            % ---- GUI mode -> graphic data is given to the appropriate handler --%
        else
            currentHandler = axesHandles.genePlot;
        end
        
        if (~isempty(dataGroups(k).colorCodes))
            
            hold(currentHandler);
            for l=1:size(dataX,2)
                % if we have more data column than color, we start over the
                % the colors
                if l < size(dataGroups(k).colorCodes,2)
                    colorCode = dataGroups(k).colorCodes{l};
                else
                    colorCode = dataGroups(k).colorCodes{mod(l,size(dataGroups(k).colorCodes,2))+1};
                end
                plot(currentHandler,t_ode/3600,dataX(:,l),colorCode);
            end
        else
            plot(currentHandler,t_ode/3600,dataX);
        end
        
        
        lgh = legend(currentHandler,listOfProteins, 'Location', 'NorthEast');
        if verLessThan('matlab', '8.4')
            legend(lgh,'boxoff');
        else
            lgh.Box = 'off';
        end
        ylabel(currentHandler,'Species amounts [nM]');
        xlabel(currentHandler,'Time [h]');
        title(currentHandler,dataGroups(k).plotName);
        
        % add the processed data to the output structure
        processedData{1} = {['Time' listOfProteins],[t_ode dataX]};
        
        %%%%%%%%%%% Resource usage plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif(strcmp(dataGroups(k).plotName,'Resource usage'))
        
        listOfResources = {'AGTP', 'CUTP', 'AA','RNAP','Ribo'};
        dataX = getDataForSpecies(currModelObj,x_ode,listOfResources);
        
        % ---- Calling the txtl_plot standalone -> figure is generated -------%
        if strcmp(operationMode,'standalone')
            currentHandler = subplot(224);
            % ---- GUI mode -> graphic data is given to the appropriate handler --%
        else
            currentHandler = axesHandles.resourceUsage;
        end
        
        mMperunit = 100 / 1000;			% convert from NTP, AA units to mM
        plot(currentHandler,...
            t_ode/3600, dataX(:, 1)/dataX(1, 1), 'b-', ...
            t_ode/3600, dataX(:, 2)/dataX(1, 2), 'r-', ...
            t_ode/3600, dataX(:, 3)/dataX(1, 3), 'b--', ...
            t_ode/3600, dataX(:, 4)/dataX(1, 4), 'r--', ...
            t_ode/3600, dataX(:, 5)/dataX(1, 5), 'g--');
        
        title(currentHandler,'Resource usage');
        lgh = legend(currentHandler,...
            {'AGTP [nM]', 'CUTP [nM]', 'AA [mM]', 'RNAP [nM]', 'Ribo [nM]'}, ...
            'Location', 'Best');
        if verLessThan('matlab', '8.4')
            legend(lgh,'boxoff');
        else
            lgh.Box = 'off';
        end
        ylabel(currentHandler,'Species amounts [normalized]');
        xlabel(currentHandler,'Time [h]');
        
        % add the processed data to the output structure
        processedData{3} = {['Time' listOfResources],[t_ode dataX]};
        
        %%%%%%%%%%% Error Handling  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        disp('no option was provided!');
        
    end % end of if dataGroups
    
end % end of for

end


function totalAmount = totalAmountOfSpecies(listOfSpecies,x_ode,SpecieString)


matchSpecie = regexp(listOfSpecies,SpecieString,'match');
if ~isempty(vertcat(matchSpecie{:}))
    totalAmount{1,1} = SpecieString;
    totalAmount{1,2} = [];
    binVec = cellfun(@(x) isempty(x),matchSpecie);
    indx = find(binVec == 0);
    totalAmount{1,2} =  sum(x_ode(:,indx),2);
    
else
    totalAmount{1,1} = SpecieString;
    totalAmount{1,2} =  0;
    warning('total concentration: no match was found for: %s',...
        regString);
    
end

end



function autoSpecies = extractRegexpAndExecute(regexp_ind,listOfSpecies,dataSource)
autoSpecies = {};
for z=1:size(regexp_ind,1)
    r_str = strrep(dataSource(z),'#','');
    specie_match = regexp(listOfSpecies,r_str,'tokens','once');
    autoSpecies{z} = vertcat(specie_match{:});
    
    %     for l=1:size(currModelObj.Species,1)
    %         specie_match = regexp(currModelObj.Species(l).Name,r_str,'tokens');
    %         if(~cellfun(@(x) isempty(x),specie_match))
    %             autoSpecies(end+1) = specie_match{1}{1};
    %         end
    %     end
end
autoSpecies = vertcat(autoSpecies{:});
end


function dataX = getDataForSpecies(currModelObj,x_ode,listOfSpecies)
% collect data for the listed species from the simulation result array (x_ode)

indexNum = findspecies(currModelObj, listOfSpecies);
notASpecie = find(indexNum == 0);
if (~isempty(notASpecie))
    for k=1:size(notASpecie,2)
        error('not valid specie name: %s!',listOfSpecies{notASpecie(k)});
    end
end
dataX = x_ode(:,indexNum);

end

function rgbVector = convertLabelToRGBValue(label)
rgbVector = rem(floor((strfind('kbgcrmyw', label) - 1) * [0.25 0.5 1]), 2);

end


function [ColorMtx,LineStyle] = getColorAndLineOrderByUserData(listOfItems)
styleOptions = regexp(listOfItems,'(.)(.*)','tokens','once');
ColorMtx = [];
LineStyle = {};

for l=1:size(styleOptions,2)
    ColorMtx(l,:) = convertLabelToRGBValue(styleOptions{1,l}{1});
    if isempty(styleOptions{1,l}{2})
        LineStyle{l} = '-';
    else
        LineStyle{l} = styleOptions{1,l}{2};
    end
    
end
end