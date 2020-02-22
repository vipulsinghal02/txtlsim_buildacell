% plotEverything
numSubplots = 6;
[~,listOfSpecies] = getstoichmatrix(Mobj{1});
numfullfigs = floor(length(listOfSpecies)/numSubplots);
numSubplotsPartialFig = mod(length(listOfSpecies),numSubplots);

cellOfSpecies = cell(numSubplots,1);


for jj = 1:numfullfigs
    for k = 1:numSubplots
        cellOfSpecies{k} = listOfSpecies{numSubplots*(jj-1)+ k};
    end
    plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, cellOfSpecies, LegendList)
end

cellOfSpecies = cell(numSubplotsPartialFig,1);
for k = 1:numSubplotsPartialFig
        cellOfSpecies{k} = listOfSpecies{numSubplots*(numfullfigs)+ k};
end
plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, cellOfSpecies, LegendList)
        
        
    