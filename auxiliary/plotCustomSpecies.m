function varargout = plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, legendList, varargin)
% plotting routine: input list of species, mobj cell array, data cell array, title containing what is being plotted and what speie is being varied, and cell array of parameter values

if class(Mobj) == 'cell'
    
    colororder1 = lines;
    colororder2 = [0.8 0 0;
        0 0.8 0;
        160/255 32/255 240/255;
        0 1 1;
        1 69/255 0;
        112/255 138/255 144/255;
        188/255 143/255 143/255
        0 0 0.8;];
    colororder3 = [colororder2;colororder1];
    numberOfSubplots = size(cellOfSpecies,1);
    
    % curves per subplot
    numberOfCurves = length(Mobj);
    
    
    outputCurves = cell(numberOfSubplots,numberOfCurves);
    notEmpty = ~cellfun(@isempty,cellOfSpecies);
    % a given curve may be a sum of other curves
    numberOfSpeciesToAdd = sum(notEmpty,2);
    figure
    for i = 1:numberOfSubplots
        subplot(numberOfSubplots, 1,i)
        
        for k = 1:numberOfCurves
            dataToPlot = 0;
            for m = 1:numberOfSpeciesToAdd(i)
                
                specieIndex = findspecies(Mobj{k}, cellOfSpecies{i,m});
                dataToPlot = dataToPlot + x_ode{k}(:,specieIndex);
            end
            h(k) = plot(t_ode{k}/60, dataToPlot);
            outputCurves{i,k} = dataToPlot;
            hold on
            set(h(k), 'Color', colororder3(k,:), 'LineWidth', 1.5);
            hold on
            
        end
%         lgh = legend(h, legendList, 'Location', 'NorthEastOutside')
%         set(lgh,'FontSize',8);
        title(titleString{i},'FontSize', 7);
    end
    
    [ax1,h1]=suplabel('time/min');
    [ax2,h2]=suplabel('conc or RFU','y');
    
    if nargin==7
        if strcmp(varargin{1},'save')
            folderdate = datestr(now,'yyyymmmmdd_HHMMSS');
            set(0,'DefaultFigureVisible','on');
            % scrsz = get(0,'ScreenSize');
            mainDir = pwd;
            saveTo = [mainDir '\examples\CSHL\' folderdate];
            mkdir(saveTo)
            %figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])%,'Name','Simulation Plot Window','NumberTitle','off')
            
            cd(saveTo)
            print('-dtiff','-r100',[titleString{1} '.tiff'])
            saveas(gcf, [titleString{1} '.fig'])
            cd(mainDir)
        end
    end

    varargout{1} = outputCurves;

else
    warning('Mobj and other things must be cell arrays')
    
    
end
