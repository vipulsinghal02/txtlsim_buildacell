function varargout = plotCustomSpecies2(mobj, x_ode, t_ode, cellofspecies, varargin)
% plotting routine: input list of species, mobj cell array, data cell array, 
% title containing what is being plotted and what speie is being varied, 
% and cell array of parameter values
numvarargs = length(varargin);
optargs = {{}, [],[],[], []};
optargs(1:numvarargs) = varargin;
[legendList, saveflag, folderName, figureName, fighandle] = optargs{:};
 scrsz = get(0,'screensize');
 if ~isempty(fighandle)
 figure(fighandle )%,'name','simulation plot window','numbertitle','off')
 set(fighandle, 'position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3]);
 
 
 else
   figure( 'position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])
   %,'name','simulation plot window','numbertitle','off')
 end
 
   
if iscell(mobj)
    
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
    
    % size of subplot grid
    numberofsubplotsy = size(cellofspecies,1);
    numberofsubplotsx = size(cellofspecies,2);
    
    % curves per subplot
    numberofcurves = length(mobj);
    
%     figure
    for i = 1:numberofsubplotsy
        for jj = 1:numberofsubplotsx
            plotnumber = (i-1)*numberofsubplotsx + jj;
            subplot(numberofsubplotsy, numberofsubplotsx,plotnumber)
            
            for k = 1:numberofcurves
                specieindex = findspecies(mobj{k}, cellofspecies{i,jj});
                datatoplot = x_ode{k}(:,specieindex);
                h(k) = plot(t_ode{k}/60, datatoplot);
                hold on
                set(h(k), 'color', colororder3(k,:), 'linewidth', 1.5);
            end
            title(cellofspecies{i,jj},'fontsize', 10);
%             
%             suplabel('time/min');
%             suplabel('conc or rfu','y');
        end
    end
    if ~isempty(legendList)
        lgh = legend(h, legendList, 'Location', 'NorthEastOutside');
        set(lgh,'FontSize',10);
        
    end
    
    if strcmp(saveflag,'save')
        
        set(0,'defaultfigurevisible','on');
        
        maindir = pwd;
        saveto = [maindir '\Sim Data\']; %folderName
        mkdir(saveto)
        
        
        cd(saveto)
        print('-dtiff','-r200',[folderName figureName '.tiff'])
        saveas(gcf, [folderName figureName '.fig'])
        cd(maindir)
    end
    
    
else
    warning('txtltoolbox:plotcustomspecies',...
        'mobj and other things must be cell arrays, nothing will be plotted')
    
    
end
