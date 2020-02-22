load('mergedExperimentFiles.mat')

colorCodes = {'r','b','g','c','m','k',[1 1 .5],[.7 .5 .2],[0 1 .2],[.35 .8 .8],[.9 0 .4],[1 .2 .2]};

RFUumConvert = 1723;   % (1723 a.u. = 1 uM)
RFUumConvertMG = 7.75;  % (7.75 a.u. = 1 nM)

disp('everything is loaded, ready to roll!');


t_vec_mins=mergedExpFile(1).t_vec/60;

%%

%%%%  GFP KINETICS
plotSelect=[4, 5, 6, 7, 8, 10];
figure('Name','GFP kinetics');
hold on
h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,11+plotSelect,2)/RFUumConvert,mergedExpFile(1).noBg_std(:,11+plotSelect,2)/RFUumConvert/sqrt(3),0.3,colorCodes);
axis([-10,850,-1000/RFUumConvert,35000/RFUumConvert])
xlabel('Time [mins]')
ylabel('deGFP (uM)')
grid on
legendStr = cellstr(num2str(mergedExpFile(1).concentrations(plotSelect)', 'pr-gfp-mg15 %0.2g nM'));
legend(legendStr,'Location','NorthWest')



%%%%%  MG KINETICS WITH CORRECT BACKGROUNDS SUBTRACTED
mgBackground = mergedExpFile(1).Data_mean(:,1:11,1);
mgNoBgMean=mergedExpFile(1).Data_mean(:,12:22,1)-mgBackground;

for i=1:11
    mgNoBgMean(:,i)=mean([expFile(1).Data(:,12+i,1)-mgBackground(:,i) expFile(2).Data(:,12+i,1)-mgBackground(:,i) expFile(2).Data(:,35+i,1)-mgBackground(:,i)],2)/RFUumConvertMG;
    mgNoBgErr(:,i)=std([expFile(1).Data(:,12+i,1)-mgBackground(:,i) expFile(2).Data(:,12+i,1)-mgBackground(:,i) expFile(2).Data(:,35+i,1)-mgBackground(:,i)],0,2)/sqrt(3)/RFUumConvertMG;
end

figure('Name','MG kinetics');
hold on
% errorbar(repmat(t_vec_mins,1,length(plotSelect)),mgNoBgMean(:,plotSelect),mgNoBgErr(:,plotSelect));
h = stdshade(t_vec_mins,mgNoBgMean(:,plotSelect),mgNoBgErr(:,plotSelect),0.3,colorCodes);
%h = stdshade(t_vec_mins,mergedExpFile(1).noBg_mean(:,plotSelect+11,1)/RFUumConvertMG,mergedExpFile(1).noBg_std(:,plotSelect+11)/RFUumConvertMG,0.3,colorCodes);
axis([-10,850,-100/RFUumConvertMG,3500/RFUumConvertMG])
xlabel('Time [mins]')
ylabel('MG fluorescence (A.U.)')
grid on
legendStr = cellstr(num2str(mergedExpFile(1).concentrations(plotSelect)', 'pr-gfp-mg15 %0.2g nM'));
legend(h,legendStr,'Location','NorthEast')