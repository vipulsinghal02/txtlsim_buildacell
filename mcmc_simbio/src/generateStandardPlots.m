function mcat_l10 = generateStandardPlots(datafiles, titlestr, legends, plotmode)
% datafiles need to be from the same batch, so that things like 
% {'simulatedDataMatrix', 'dosedInitVals','measuredSpecies',...
%     'exportedMdlObj','tspan'};
% are compatible


if strcmp(plotmode, 'txtl')
    vars = {'simulatedDataMatrix', 'dosedInitVals','measuredSpecies',...
    'exportedMdlObj','tspan'};
    load(datafiles{1}, vars{:});
    mcat = catMC(datafiles);
    mcat_l10 = mcat/log(10);
    plotEstimTraces003(mcat,exportedMdlObj,tspan, ...
    simulatedDataMatrix, dosedInitVals,...
    measuredSpecies);% , 'paramID', [1 2 3]
elseif strcmp(plotmode, 'simplemodel')
    
    mcat = catMC(datafiles);
    mcat_l10 = mcat/log(10);
end

clear mcat

plotChains(mcat_l10, 100, legends);
suptitle(titlestr)

figure
[C,lags,ESS]=eacorr(mcat_l10);
plot(lags,C,'.-',lags([1 end]),[0 0],'k');
grid on
xlabel('lags')
ylabel('autocorrelation');
text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),...
    'verticalalignment','bottom','horizontalalignment','right')
suptitle(titlestr)
% 
figure;
ecornerplot_vse(mcat_l10,...
    'ess', 30,'ks',true, 'color',[.6 .35 .3], ...
    'names', legends, 'fontsize', 16)
% suptitle(titlestr)

figure;
ecornerplot_vse(mcat_l10,...
    'scatter', true,'transparency',0.01, 'color',[.6 .35 .3], ...
    'names', legends);
suptitle(titlestr);





end

