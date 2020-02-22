%%
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

function dataGroup = txtl_getDefaultPlotDataStruct()

dataGroup(1).plotName = 'DNA and mRNA';
dataGroup(1).SpeciesToPlot = {'ALL_DNA'};
dataGroup(1).colorCodes = {'b','r','g','b--','r--','g--','c','y','w','k'}; 

% Gene Expression Plot
dataGroup(2).plotName = 'Gene Expression';
dataGroup(2).SpeciesToPlot = {'ALL_PROTEIN'};
dataGroup(2).colorCodes = {'b','r','g','b--','r--','g--','c','y','w','k'}; 

% Resource Plot
dataGroup(3).plotName = 'Resource usage';

end
