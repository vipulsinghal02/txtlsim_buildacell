clear variables;
close all
workingDir = 'mcmc_simbio/exp_data/public_data/';

% run merger
%all_data_merger;
load([workingDir 'mergedExperimentFiles.mat'])

colorCodes = {'r','b','g','c','m','k',[1 1 .5],[.7 .5 .2],[0 1 .2],[.35 .8 .8],[.9 0 .4],[1 .2 .2]};

RFUumConvert = 1723;   % (1723 a.u. = 1 uM)
RFUumConvertMG = 7.75;  % (7.75 a.u. = 1 nM)

disp('everything is loaded, ready to roll!');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr-gfp vs Pr-gfp-mg15

prgfp_prgfpmg15_dataMtx = reshape(mergedExpFile(1).noBg_mean(end,:,2),11,2);
prgfp_prgfpmg15_stdMtx  = reshape(mergedExpFile(1).noBg_std(end,:,2),11,2);

t_vec_mins=mergedExpFile(1).t_vec/60;

%%%%  BAR PLOTS COMPARING GFP WITH GFP-MG
% small_amounts = 1:5;
% large_amounts = 6:11;
% figure('Name','endpoint small')
% barweb(prgfp_prgfpmg15_dataMtx(small_amounts,:),prgfp_prgfpmg15_stdMtx(small_amounts,:),[],mergedExpFile(1).concentrations(small_amounts),'GFP protein expression endpoints ',[],[],'summer','y',{'pr-gfp','pr-gfp-mg15'},[],'plot');
% xlabel('plasmid DNA concentration [nM]')
% ylabel('GFP r.f.u.')
% figure('Name','endpoint large')
% barweb(prgfp_prgfpmg15_dataMtx(large_amounts,:),prgfp_prgfpmg15_stdMtx(large_amounts,:),[],mergedExpFile(1).concentrations(large_amounts),'GFP protein expression endpoints ',[],[],'summer','y',{'pr-gfp','pr-gfp-mg15'},[],'plot');
% xlabel('plasmid DNA concentration [nM]')
% ylabel('GFP r.f.u.')
figure('Name','barplot log scale')
barweb(prgfp_prgfpmg15_dataMtx(2:end-1,:),prgfp_prgfpmg15_stdMtx(2:end-1,:),...
    [],mergedExpFile(1).concentrations(2:end-1),'GFP protein expression endpoints ',...
    [],[],'summer','y',{'pr-gfp','pr-gfp-mg15'},[],'plot');

ylim([[50 50000]])
h=gca;
set(h,'YScale','log');

%%%%  GFP ENDPOINT VS CONCENTRATION
% linear fit
f = fittype('poly1');
selectedConc = 1:6;
x = mergedExpFile(1).concentrations([selectedConc]+1)';
y = mergedExpFile(1).noBg_mean(end,[selectedConc]+12,2)'/RFUumConvert;
fitobject = fit(x,y,f);
figure('Name','endpoint vs concentration')
hold on
errorbar(mergedExpFile(1).concentrations(2:end),mergedExpFile(1).noBg_mean(end,13:22,2)/RFUumConvert,mergedExpFile(1).noBg_std(end,13:22,2)/RFUumConvert/sqrt(3))
% h = stdshade(mergedExpFile(1).concentrations(2:end)',mergedExpFile(1).noBg_mean(end,13:22,2)',(mergedExpFile(1).noBg_std(end,13:22,2)/sqrt(3))',0.3,'b');
plot(logspace(-1,0.8),fitobject(logspace(-1,0.8)),'r')
axis([0.03,50,-1000/RFUumConvert,35000/RFUumConvert])
xlabel('plasmid DNA concentration (nM)')
ylabel('Endpoint deGFP (uM)')
set(gca,'XScale','log','YScale','linear')
%grid on
hold off
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
legend(h,legendStr,'Location','NorthWest')
%%
%%%%%  MG V1.0
% plot(mergedExpFile(1).t_vec/60,mergedExpFile(1).noBg_mean(:,[16:end],1),'LineWidth',1.7)
% xlabel('Time [min]')
% ylabel('MGApt r.f.u.')
% title('E15 T=29C Ex/Em: 610/655 Gain: 150, MGDye: ~10uM');
% legend('pr-gfp-mg 0.5nM','pr-gfp-mg 1nM','pr-gfp-mg 2nM','pr-gfp-mg 5nM','pr-gfp-mg 10nM','pr-gfp-mg 20nM','pr-gfp-mg 31.5nM')
% yL = get(gca,'YLim');
%  line([200 200],yL,'Color','c');
%  text(200,800,...
%       '\leftarrow t= 200 min',...
%       'FontSize',12)

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
%%
%%%%%  TOTAL MG
totalMGmean=sum(mergedExpFile(1).noBg_mean(:,12:end,1));%mean([sum(expFile(1).Data(:,13:23,1)-mgBackground); sum(expFile(2).Data(:,13:23,1)-mgBackground); sum(expFile(2).Data(:,36:46,1)-mgBackground)]);
totalMGerr= sum(mergedExpFile(1).noBg_std(:,12:end,1));%std([sum(expFile(1).Data(:,13:23,1)-mgBackground); sum(expFile(2).Data(:,13:23,1)-mgBackground); sum(expFile(2).Data(:,36:46,1)-mgBackground)])/sqrt(3);

f = fittype('poly1');
selectedConc = 1:5;
x = mergedExpFile(1).concentrations(selectedConc)';
y = totalMGmean(selectedConc)';
fitobject = fit(x,y,f);
figure('Name','total RNA produced vs concentration')
hold on
errorbar(mergedExpFile(1).concentrations(2:end),totalMGmean(2:end),totalMGerr(2:end))
%errorbar(mergedExpFile(1).concentrations(2:end),totalMGmean(2:end),totalMGerr(2:end))
%plot(mergedExpFile(1).concentrations(2:end),fitobject(mergedExpFile(1).concentrations(2:end)),'r')
axis([0.03,50,-15000,400000])
xlabel('plasmid DNA concentration (nM)')
ylabel('Integrated MG (A.U.)')
set(gca,'XScale','log','YScale','linear')
%grid on
hold off
%%
%%%%%  TOTAL MG VS GFP
GFPendpoints=mergedExpFile(1).noBg_mean(end,12:22,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(1).noBg_std(end,12:22,2)/RFUumConvert/sqrt(3);
figure;
hold on
% h=plot(totalMGmean,GFPendpoints,'o');
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
hverrorbar(totalMGmean(1:11),GFPendpoints,totalMGerr(1:11),GFPendpoints_err,'b')

ylim([[-1000/RFUumConvert 35000/RFUumConvert]])
xlim([[-10000 490000]])
grid on
hold off
ylabel('GFP r.f.u.')
xlabel('total MG r.f.u.')



e15_specs_0902 = importdata('ACS_paper_data/130902MG.txt');

e15_specs_0902.noBG = e15_specs_0902.data(2:2:end,:) - repmat(e15_specs_0902.data(2,:),size(e15_specs_0902.data(2:2:end,:),1),1);
e15_specs_0902.noBG = e15_specs_0902.noBG -repmat(e15_specs_0902.noBG(:,1),1,size(e15_specs_0902.noBG(2:2:end,:),2));
e15_specs_0902.normalized = e15_specs_0902.noBG./repmat(max(e15_specs_0902.noBG')',1,size(e15_specs_0902.noBG,2));
e15_specs_0902.t_vec = e15_specs_0902.data(1,:);

ll = e15_specs_0902.noBG(2:end,:)./repmat(max(e15_specs_0902.noBG(2:end,:)')',1,size(e15_specs_0902.noBG(2:end,:),2))
kk = mgNoBgMean(:,[6:7])./repmat(max(mgNoBgMean(:,[6:7])),size(mgNoBgMean(:,[6:7]),1),1);
figure('Name','Compare spex vs biotek')
hold on
plot(e15_specs_0902.t_vec/60,ll,'--');
plot(t_vec_mins,kk);
% h = stdshade(t_vec_mins,mgNoBgMean(:,[6:7]),mgNoBgErr(:,[6:7]),0.3,colorCodes);
hold off
legend('pr-gfp-mg15 1nM spex','pr-gfp-mg15 2nM spex','pr-gfp-mg15 3nM spex','pr-gfp-mg15 1nM biotek','pr-gfp-mg15 2nM biotek')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr-gfp-mg15 + NTP

prgfpmg15_NTP_DataMtx = reshape(mergedExpFile(2).noBg_mean(end,:,2),10,1);
prgfpmg15_NTP_StdMtx  = reshape(mergedExpFile(2).noBg_std(end,:,2),10,1);
t_vec_mins=mergedExpFile(2).t_vec/60;

plotSelect=[4, 5, 6, 7, 8, 10];

% small_amounts = 1:5;
% large_amounts = 6:10;
% figure('Name','endpoint small')
% barweb(prgfpmg15_NTP_DataMtx(small_amounts,:),prgfpmg15_NTP_StdMtx(small_amounts,:),[],mergedExpFile(3).constructNames,'GFP protein expression endpoints ',[],[],'summer','y',mergedExpFile(3).concentrations(small_amounts),[],'plot');
% figure('Name','endpoint large')
% barweb(prgfpmg15_NTP_DataMtx(large_amounts,:),prgfpmg15_NTP_StdMtx(large_amounts,:),[],mergedExpFile(3).constructNames,'GFP protein expression endpoints ',[],[],'summer','y',mergedExpFile(3).concentrations(large_amounts),[],'plot');
%
% plot(mergedExpFile(2).t_vec/60,mergedExpFile(2).noBg_mean(:,5:end,1),'LineWidth',1.7);
% xlabel('Time [min]')
% ylabel('MGApt r.f.u.')
% title('E15 T=29C Ex/Em: 610/655 Gain: 150, MGDye: ~10uM');
% legend('pr-gfp-mg 0.5nM +NTP','pr-gfp-mg 1nM +NTP','pr-gfp-mg 2nM +NTP','pr-gfp-mg 5nM +NTP','pr-gfp-mg 10nM +NTP','pr-gfp-mg 20nM +NTP')
% yL = get(gca,'YLim');
%  line([200 200],yL,'Color','c');
%  text(200,800,...
%       '\leftarrow t= 200 min',...
%       'FontSize',12)

%%%%  GFP+NTP KINETICS
figure('Name','GFP kinetics');
hold on
h = stdshade(t_vec_mins,mergedExpFile(2).noBg_mean(:,plotSelect,2)/RFUumConvert,mergedExpFile(2).noBg_std(:,plotSelect,2)/RFUumConvert/sqrt(3),0.3,colorCodes);
axis([-10,850,-1000/RFUumConvert,35000/RFUumConvert])
xlabel('Time (mins)')
ylabel('deGFP fluorescence (A.U.)')
grid on
legendStr = cellstr(num2str(mergedExpFile(2).concentrations(plotSelect)', 'Pr-deGFP-MGapt %0.2g nM'));
legend(h,legendStr,'Location','NorthWest')

%%%%  GFP+NTP ENDPOINT VS CONCENTRATION
% linear fit
f = fittype('poly1');
selectedConc = 1:7;
x = mergedExpFile(2).concentrations([selectedConc])';
y = mergedExpFile(2).noBg_mean(end,[selectedConc],2)';
fitobject = fit(x,y,f);
figure('Name','endpoint vs concentration')
hold on
h = errorbar(mergedExpFile(2).concentrations,mergedExpFile(2).noBg_mean(end,:,2),mergedExpFile(2).noBg_std(end,:,2)/sqrt(3));
% h = stdshade(mergedExpFile(3).concentrations',mergedExpFile(2).noBg_mean(end,:,2)',(mergedExpFile(2).noBg_std(end,:,2)/sqrt(3))',0.3,'b');
plot(0:0.2:2.8,fitobject(0:0.2:2.8),'r')
axis([-1,23,-1000,27000])
xlabel('plasmid DNA concentration [nM]')
ylabel('GFP r.f.u.')
%set(gca,'XScale','log','YScale','linear')
grid on
hold off


%%%%% MG + NTP WITH CORRECT BACKGROUNDS SUBTRACTED

t_vec_mins=mergedExpFile(2).t_vec/60;

% mgNTPNoBg=mergedExpFile(2).Data_mean(:,:,1)-mergedExpFile(3).Data_mean(:,:,1);
% figure;
% plot(t_vec_mins,mgNTPNoBg,'LineWidth',1.7);

% mgNTPBackground = mergedExpFile(3).Data_mean(:,:,1);
% for i=1:10
%     mgNTPNoBgMean(:,i)=mean([expFile(3).Data(:,1+i,1)-mgNTPBackground(:,i) expFile(3).Data(:,11+i,1)-mgNTPBackground(:,i) expFile(3).Data(:,21+i,1)-mgNTPBackground(:,i)],2);
%     mgNTPNoBgErr(:,i)=std([expFile(3).Data(:,1+i,1)-mgNTPBackground(:,i) expFile(3).Data(:,11+i,1)-mgNTPBackground(:,i) expFile(3).Data(:,21+i,1)-mgNTPBackground(:,i)],0,2)/sqrt(3);
% end
figure('Name','MG kinetics');
hold on
% errorbar(repmat(t_vec_mins,1,length(plotSelect)),mgNoBgMean(:,plotSelect),mgNoBgErr(:,plotSelect));
% h = stdshade(t_vec_mins,mgNTPNoBgMean(:,plotSelect),mgNTPNoBgErr(:,plotSelect),0.3,colorCodes);
h = stdshade(t_vec_mins,mergedExpFile(2).noBg_mean(:,plotSelect),mergedExpFile(2).noBg_std(:,plotSelect)/sqrt(3),0.3,colorCodes);
axis([-10,850,-100,3500])
xlabel('Time (mins)')
ylabel('MGapt fluorescence (A.U.)')
grid on
legendStr = cellstr(num2str(mergedExpFile(1).concentrations(plotSelect)', 'Pr-deGFP-MGapt %0.2g nM'));
legend(h,legendStr,'Location','NorthEast')

%%%%%  TOTAL MG VS GFP

% totalMGmean=mean([sum(expFile(3).Data(:,2:11,1)-mgNTPBackground); sum(expFile(3).Data(:,12:21,1)-mgNTPBackground); sum(expFile(3).Data(:,22:31,1)-mgNTPBackground)]);
% totalMGerr=std([sum(expFile(3).Data(:,2:11,1)-mgNTPBackground); sum(expFile(3).Data(:,12:21,1)-mgNTPBackground); sum(expFile(3).Data(:,22:31,1)-mgNTPBackground)])/sqrt(3);
totalMGmean=sum(mergedExpFile(2).noBg_mean(:,:,1));
totalMGerr=sum(mergedExpFile(2).noBg_std(:,:,1)/sqrt(3));
% ZT: NOTE that std is different, because the order of operations!

GFPendpoints=mergedExpFile(2).noBg_mean(end,:,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(2).noBg_std(end,:,2)/RFUumConvert/sqrt(3);
figure;
hold on
% h=plot(totalMGmean,GFPendpoints,'o');
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
hverrorbar(totalMGmean,GFPendpoints,totalMGerr,GFPendpoints_err,'b')

ylim([[-1000/RFUumConvert 35000/RFUumConvert]])
xlim([[-10000 490000]])
grid on
hold off
ylabel('GFP r.f.u.')
xlabel('total MG r.f.u.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr1-gfp-mg15

pr1DataMtx = reshape(mergedExpFile(11).noBg_mean(end,:,2),10,1);
pr1StdMtx  = reshape(mergedExpFile(11).noBg_std(end,:,2),10,1);

t_vec_mins=mergedExpFile(11).t_vec/60;
plotSelect=[4, 5, 6, 7, 8, 10]-1;

small_amounts = 1:5;
large_amounts = 6:10;
figure('Name','endpoint small')
legendStr = cellstr(num2str(mergedExpFile(11).concentrations', 'Pr1-deGFP-MGapt %0.2g nM'));

barweb(pr1DataMtx(small_amounts,:),pr1StdMtx(small_amounts,:),[],mergedExpFile(11).constructNames,'GFP protein expression endpoints ',[],[],'summer','y',legendStr(small_amounts),[],'plot');
figure('Name','endpoint large')
barweb(pr1DataMtx(large_amounts,:),pr1StdMtx(large_amounts,:),[],mergedExpFile(11).constructNames,'GFP protein expression endpoints ',[],[],'summer','y',legendStr(large_amounts),[],'plot');

%%%%%  Pr1 GFP kinetics
figure('Name','Pr1 GFP kinetics');
hold on
h = stdshade(t_vec_mins,mergedExpFile(11).noBg_mean(:,1:10,2)/RFUumConvert,mergedExpFile(11).noBg_std(:,1:10,2)/RFUumConvert/sqrt(3),0.3,colorCodes);
axis([-10,850,-1000/RFUumConvert,35000/RFUumConvert])
xlabel('Time (mins)')
ylabel('deGFP fluorescence (A.U.)')
grid on
legendStr = cellstr(num2str(mergedExpFile(11).concentrations(1:10)', 'Pr1-deGFP-MGapt %0.2g nM'));
legend(h,legendStr,'Location','Best')

%%%%%  TOTAL MG VS GFP
%use the background from Pr-GFP case:
mgBackground = mergedExpFile(1).Data_mean(:,1:10,1);

totalMGmean=mean([sum(expFile(4).Data(:,2:11,1)-mgBackground); sum(expFile(4).Data(:,12:21,1)-mgBackground); sum(expFile(4).Data(:,22:31,1)-mgBackground)]);
totalMGerr=std([sum(expFile(4).Data(:,2:11,1)-mgBackground); sum(expFile(4).Data(:,12:21,1)-mgBackground); sum(expFile(4).Data(:,22:31,1)-mgBackground)])/sqrt(3);

GFPendpoints=mergedExpFile(11).noBg_mean(end,:,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(11).noBg_std(end,:,2)/RFUumConvert/sqrt(3);
figure;
hold on
% h=plot(totalMGmean,GFPendpoints,'o');
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
hverrorbar(totalMGmean,GFPendpoints,totalMGerr,GFPendpoints_err,'b')
ylim([[-1000/RFUumConvert 35000/RFUumConvert]])
xlim([[-10000 490000]])
grid on
hold off
ylabel('GFP r.f.u.')
xlabel('total MG r.f.u.')

%%%%  GFP ENDPOINT VS CONCENTRATION
% linear fit
f = fittype('poly1');
x = mergedExpFile(11).concentrations(4:end-1)';
y = mergedExpFile(11).noBg_mean(end,4:end-1,2)'/RFUumConvert;
fitobject = fit(x,y,f);
figure('Name','endpoint vs concentration')
hold on
errorbar(mergedExpFile(11).concentrations(2:end),mergedExpFile(11).noBg_mean(end,2:end,2)/RFUumConvert,mergedExpFile(11).noBg_std(end,2:end,2)/RFUumConvert/sqrt(3))
% h = stdshade(mergedExpFile(1).concentrations(2:end)',mergedExpFile(1).noBg_mean(end,13:22,2)',(mergedExpFile(1).noBg_std(end,13:22,2)/sqrt(3))',0.3,'b');
plot(logspace(-1,1.2),fitobject(logspace(-1,1.2)),'r')
axis([0.03,50,-.1,6.1])
xlabel('plasmid DNA concentration (nM)')
ylabel('Endpoint deGFP (uM)')
set(gca,'XScale','log','YScale','linear')
%grid on
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr2-gfp-mg15

pr2DataMtx = reshape(mergedExpFile(4).noBg_mean(end,:,2),8,1);
pr2StdMtx  = reshape(mergedExpFile(4).noBg_std(end,:,2),8,1);

t_vec_mins=mergedExpFile(4).t_vec/60;
plotSelect=[4, 5, 6, 7, 8, 10]-2;

small_amounts = 1:4;
large_amounts = 5:8;
figure('Name','endpoint small')
barweb(pr2DataMtx(small_amounts,:),pr2StdMtx(small_amounts,:),[],mergedExpFile(4).constructNames,'GFP protein expression endpoints ',[],[],'summer','y',mergedExpFile(4).concentrations(small_amounts),[],'plot');
figure('Name','endpoint large')
barweb(pr2DataMtx(large_amounts,:),pr2StdMtx(large_amounts,:),[],mergedExpFile(4).constructNames,'GFP protein expression endpoints ',[],[],'summer','y',mergedExpFile(4).concentrations(large_amounts),[],'plot');

%%%%%  Pr2 GFP kinetics
figure('Name','Pr2 GFP kinetics');
hold on
h = stdshade(t_vec_mins,mergedExpFile(4).noBg_mean(:,plotSelect,2)/RFUumConvert,mergedExpFile(4).noBg_std(:,plotSelect,2)/RFUumConvert/sqrt(3),0.3,colorCodes);
axis([-10,850,-1000/RFUumConvert,35000/RFUumConvert/10])
xlabel('Time (mins)')
ylabel('deGFP fluorescence (A.U.)')
grid on
legendStr = cellstr(num2str(mergedExpFile(4).concentrations(plotSelect)', 'Pr2-deGFP-MGapt %0.2g nM'));
legend(h,legendStr,'Location','NorthWest')


%%%%%  TOTAL MG VS GFP
%use the background from Pr-GFP case:
mgBackground = mergedExpFile(1).Data_mean(:,3:10,1);

totalMGmean=mean([sum(expFile(5).Data(:,2:9,1)-mgBackground); sum(expFile(5).Data(:,10:17,1)-mgBackground); sum(expFile(5).Data(:,18:25,1)-mgBackground)]);
totalMGerr=std([sum(expFile(5).Data(:,2:9,1)-mgBackground); sum(expFile(5).Data(:,10:17,1)-mgBackground); sum(expFile(5).Data(:,18:25,1)-mgBackground)])/sqrt(3);

GFPendpoints=mergedExpFile(4).noBg_mean(end,:,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(4).noBg_std(end,:,2)/RFUumConvert/sqrt(3);
figure;
hold on
% h=plot(totalMGmean,GFPendpoints,'o');
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
hverrorbar(totalMGmean,GFPendpoints,totalMGerr,GFPendpoints_err,'b')
ylim([[-.05 .5]])
xlim([[-10000 100000]])
grid on
hold off
ylabel('GFP r.f.u.')
xlabel('total MG r.f.u.')


% linear fit
f = fittype('poly1');
x = mergedExpFile(4).concentrations(4:end-1)';
y = mergedExpFile(4).noBg_mean(end,4:end-1,2)'/RFUumConvert;
fitobject = fit(x,y,f);
figure('Name','endpoint vs concentration')
hold on
errorbar(mergedExpFile(4).concentrations(2:end),mergedExpFile(4).noBg_mean(end,2:end,2)/RFUumConvert,mergedExpFile(4).noBg_std(end,2:end,2)/RFUumConvert/sqrt(3))
% h = stdshade(mergedExpFile(1).concentrations(2:end)',mergedExpFile(1).noBg_mean(end,13:22,2)',(mergedExpFile(1).noBg_std(end,13:22,2)/sqrt(3))',0.3,'b');
plot(logspace(-1,1.35),fitobject(logspace(-1,1.35)),'r')
axis([0.03,50,-.015,0.375])
xlabel('plasmid DNA concentration (nM)')
ylabel('Endpoint deGFP (uM)')
set(gca,'XScale','log','YScale','linear')
%grid on
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr2,Pr1,Pr - T7RNAP -- T7-gfp-mg15

% rudimentary barplot output for basic checks
% t7DataMtx = reshape(mergedExpFile(5).noBg_mean(end,:,2),12,3)';
% t7StdMtx  = reshape(mergedExpFile(5).noBg_std(end,:,2),12,3)';
%
% prt7 = reshape(t7DataMtx(1,:),3,4)';
% prt7_std = reshape(t7StdMtx(1,:),3,4)';
% figure('Name','pr')
% legendStr = cellstr(num2str(mergedExpFile(5).concentrations{1}', '%0.2g nM'));
% legendStr2 = cellstr(num2str(mergedExpFile(5).concentrations{2}', '%0.2g nM'));
% barweb(prt7,prt7_std,[],legendStr,'GFP protein expression endpoints ',[],[],'summer','y',legendStr2,[],'plot');
%
% pr1t7 = reshape(t7DataMtx(2,:),3,4)';
% pr1t7_std = reshape(t7StdMtx(2,:),3,4)';
% figure('Name','pr1')
% barweb(pr1t7,pr1t7_std,[],legendStr,'GFP protein expression endpoints ',[],[],'summer','y',legendStr2,[],'plot');
%
% pr2t7 = reshape(t7DataMtx(3,:),3,4)';
% pr2t7_std = reshape(t7StdMtx(3,:),3,4)';
% figure('Name','pr2')
% barweb(pr2t7,pr2t7_std,[],legendStr,'GFP protein expression endpoints ',[],[],'summer','y',legendStr2,[],'plot');

%%%%%  MG KINETICS WITH CORRECT BACKGROUNDS SUBTRACTED
% mgBackground = repmat(mergedExpFile(1).Data_mean(:,[6 7 9],1),1,4);
t_vec_mins=mergedExpFile(5).t_vec(1:233)/60;
timeRange=[1:233];
plotSelect=[1, 2, 3];

% for i=1:12
%     mgNoBgMean(:,i)=mean([expFile(6).Data(:,2+i,1)-mgBackground(:,i) expFile(7).Data(:,2+i,1)-mgBackground(:,i) expFile(8).Data(:,2+i,1)-mgBackground(:,i)],2)/RFUumConvertMG;
%     mgNoBgErr(:,i)=std([expFile(6).Data(:,2+i,1)-mgBackground(:,i) expFile(7).Data(:,2+i,1)-mgBackground(:,i) expFile(8).Data(:,2+i,1)-mgBackground(:,i)],0,2)/sqrt(3)/RFUumConvertMG;
% end

for k=0:3
    figure('Name',sprintf('MG kinetics, %0.2g nM Pr-T7RNAP',mergedExpFile(5).concentrations{1}(k+1)));
    hold on
    idx = plotSelect+k*3;
%      h = stdshade(t_vec_mins,mgNoBgMean(1:233,idx),mgNoBgErr(1:233,plotSelect),0.3,colorCodes);
 h = stdshade(t_vec_mins,mergedExpFile(5).noBg_mean(timeRange,idx,1)/RFUumConvertMG,mergedExpFile(5).noBg_std(timeRange,idx,1)/sqrt(3)/RFUumConvertMG,0.3,colorCodes);
    xlabel('Time [mins]')
    ylabel('MGapt (nM)')
    legendStr = cellstr(num2str(mergedExpFile(5).concentrations{2}', 'PT7-deGFP-MGapt %0.2g nM'));
    legend(h,legendStr,'Location','NorthEast')
    ylim([-20 400]);
    
    
    %%%%%  GFP KINETICS WITH CORRECT BACKGROUNDS SUBTRACTED
    figure('Name',sprintf('deGFP kinetics, %0.2g nM Pr-T7RNAP',mergedExpFile(5).concentrations{1}(k+1)));
    hold on
    h = stdshade(t_vec_mins,mergedExpFile(5).noBg_mean(timeRange,idx,2)/RFUumConvert,mergedExpFile(5).noBg_std(timeRange,idx,2)/RFUumConvert/sqrt(3),0.3,colorCodes);
    xlabel('Time [mins]')
    ylabel('deGFP (uM)')
    legendStr = cellstr(num2str(mergedExpFile(5).concentrations{2}', 'PT7-deGFP-MGapt %0.2g nM'));
    legend(h,legendStr,'Location','NorthEast')
    ylim([-.5 11]);
    
end

%%%%%  TOTAL MG VS GFP
%just use Pr-GFP as background for now:
mgBackground = repmat(mergedExpFile(1).Data_mean(:,[6 7 9],1),1,4);
t_vec_mins=mergedExpFile(5).t_vec(1:233)/60;
totalTime=t_vec_mins(end);

%%%%%% Pr
totalMGmean=mean([sum(expFile(6).Data(:,3:14,1)-mgBackground); sum(expFile(7).Data(:,3:14,1)-mgBackground); sum(expFile(8).Data(:,3:14,1)-mgBackground)])*3/RFUumConvertMG/60;
totalMGerr=std([sum(expFile(6).Data(:,3:14,1)-mgBackground); sum(expFile(7).Data(:,3:14,1)-mgBackground); sum(expFile(8).Data(:,3:14,1)-mgBackground)])/sqrt(3)*3/RFUumConvertMG/60;

GFPendpoints=mergedExpFile(5).noBg_mean(end,1:12,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(5).noBg_std(end,1:12,2)/RFUumConvert/sqrt(3);
figure;
% ylim([[-1000/RFUumConvert 20000/RFUumConvert]])
% xlim([[-10000 340000]])
hold on
% h=plot(avgMGmean,GFPendpoints,'o');
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')

for k=0:3
    idx = [1:3]+k*3;
    h(k+1) = hverrorbar(totalMGmean(idx),GFPendpoints(idx),totalMGerr(idx),GFPendpoints_err(idx),colorCodes{k+1});
end

ylabel('deGFP (uM)')
xlabel('total MG r.f.u.')
grid on
ylim([[-.5 12]])
xlim([[-10000 350000]]*3/RFUumConvertMG/60)

%%%%%% Pr1
totalMGmean=mean([sum(expFile(6).Data(:,15:26,1)-mgBackground); sum(expFile(7).Data(:,15:26,1)-mgBackground); sum(expFile(8).Data(:,15:26,1)-mgBackground)])*3/RFUumConvertMG/60;
totalMGerr=std([sum(expFile(6).Data(:,15:26,1)-mgBackground); sum(expFile(7).Data(:,15:26,1)-mgBackground); sum(expFile(8).Data(:,15:26,1)-mgBackground)])/sqrt(3)*3/RFUumConvertMG/60;

GFPendpoints=mergedExpFile(5).noBg_mean(end,13:24,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(5).noBg_std(end,13:24,2)/RFUumConvert/sqrt(3);

figure;
hold on
%ylim([[-1000/RFUumConvert 20000/RFUumConvert]])
%xlim([[-10000 340000]])
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
for k=0:3
    idx = [1:3]+k*3;
    h(k+1) = hverrorbar(totalMGmean(idx),GFPendpoints(idx),totalMGerr(idx),GFPendpoints_err(idx),colorCodes{k+1});
end
hold off
ylabel('deGFP (uM)')
xlabel('total MG r.f.u.')
grid on
ylim([[-.5 12]])
xlim([[-10000 350000]]*3/RFUumConvertMG/60)

%%%%%% Pr2
totalMGmean=mean([sum(expFile(6).Data(:,27:38,1)-mgBackground); sum(expFile(7).Data(:,27:38,1)-mgBackground); sum(expFile(8).Data(:,27:38,1)-mgBackground)])*3/RFUumConvertMG/60;
totalMGerr=std([sum(expFile(6).Data(:,27:38,1)-mgBackground); sum(expFile(7).Data(:,27:38,1)-mgBackground); sum(expFile(8).Data(:,27:38,1)-mgBackground)])/sqrt(3)*3/RFUumConvertMG/60;

GFPendpoints=mergedExpFile(5).noBg_mean(end,25:36,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(5).noBg_std(end,25:36,2)/RFUumConvert/sqrt(3);

figure;
hold on
%ylim([[-1000/RFUumConvert 20000/RFUumConvert]])
%xlim([[-10000 340000]])
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
for k=0:3
    idx = [1:3]+k*3;
    h(k+1) = hverrorbar(totalMGmean(idx),GFPendpoints(idx),totalMGerr(idx),GFPendpoints_err(idx),colorCodes{k+1});
end
hold off
ylabel('deGFP (uM)')
xlabel('total MG r.f.u.')
grid on
ylim([[-.5 12]])
xlim([[-10000 350000]]*3/RFUumConvertMG/60)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr2,Pr1,Pr - T7RNAP -- T7-gfp-mg15 NTP

t7DataMtx = reshape(mergedExpFile(6).noBg_mean(end,:,2),3,4)';
t7StdMtx  = reshape(mergedExpFile(6).noBg_std(end,:,2),3,4)';

% figure('Name','pr-t7rnap + NTP')
% legendStr = cellstr(num2str(mergedExpFile(6).concentrations{1}', '%0.2g nM'));
% legendStr2 = cellstr(num2str(mergedExpFile(6).concentrations{2}', '%0.2g nM'));
% barweb(t7DataMtx,t7StdMtx,[],legendStr,'GFP protein expression endpoints ',[],[],'summer','y',legendStr2,[],'plot');


%%%%%  INTEGRATED MG VS GFP
%the background is used from Pr-GFP + NTP case:
totalTime=mergedExpFile(6).t_vec(end)/60;

totalMGmean=sum(mergedExpFile(6).noBg_mean(:,:,1))*3/RFUumConvertMG/60;
totalMGerr=sum(mergedExpFile(6).noBg_std(:,:,1)/sqrt(3))*3/RFUumConvertMG/60;

GFPendpoints=mergedExpFile(6).noBg_mean(end,:,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(6).noBg_std(end,:,2)/RFUumConvert/sqrt(3);

%%% Pr
figure;
hold on
for k=0:3
    idx = [1:3]+k*3;
    h(k+1) = hverrorbar(totalMGmean(idx),GFPendpoints(idx),totalMGerr(idx),GFPendpoints_err(idx),colorCodes{k+1});
end

grid on
hold off
ylabel('GFP uM')
xlabel('integrated MGapt nM.hr')
legend(h,'Pr-T7RNAP 0.1 nM + NTP','Pr-T7RNAP 0.2 nM + NTP','Pr-T7RNAP 0.5 nM + NTP','Pr-T7RNAP 1 nM + NTP')

%%%%%  KINETICS

t_vec_mins=mergedExpFile(6).t_vec/60;
plotSelect=[1, 2, 3];

for k=0:3
    figure('Name',sprintf('MG kinetics, %0.2g nM Pr-T7RNAP + NTP',mergedExpFile(6).concentrations{1}(k+1)));
    hold on
    idx = plotSelect+k*3;
%   h = stdshade(t_vec_mins,mgNoBgMean(1:233,idx),mgNoBgErr(1:233,plotSelect),0.3,colorCodes);
    h = stdshade(mergedExpFile(6).t_vec/60,mergedExpFile(6).noBg_mean(:,idx,1)/RFUumConvertMG,mergedExpFile(6).noBg_std(:,idx,1)/sqrt(3)/RFUumConvertMG,0.3,colorCodes);
    xlabel('Time (mins)')
    ylabel('MGapt (nM)')
    legendStr = cellstr(num2str(mergedExpFile(6).concentrations{2}', 'PT7-deGFP-MGapt %0.2g nM + NTP'));
    legend(h,legendStr,'Location','NorthEast')
    ylim([-20 900]);
end

for k=0:3
    figure('Name',sprintf('deGFP kinetics, %0.2g nM Pr-T7RNAP + NTP',mergedExpFile(6).concentrations{1}(k+1)));
    hold on
    idx = plotSelect+k*3;
    h = stdshade(t_vec_mins,mergedExpFile(6).noBg_mean(:,idx,2)/RFUumConvert,mergedExpFile(6).noBg_std(:,idx,2)/RFUumConvert/sqrt(3),0.3,colorCodes);
    xlabel('Time (mins)')
    ylabel('deGFP (uM)')
    legendStr = cellstr(num2str(mergedExpFile(6).concentrations{2}', 'PT7-deGFP-MGapt %0.2g nM + NTP'));
    legend(h,legendStr,'Location','NorthEast')
    ylim([-.5 16]);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr-s70 pr-GFP-mg15

s70DataMtx = reshape(mergedExpFile(8).noBg_mean(end,:,2),3,4)';
s70StdMtx  = reshape(mergedExpFile(8).noBg_std(end,:,2),3,4)';

figure('Name','endpoint GFP')
legendStr = cellstr(num2str(mergedExpFile(8).concentrations{1}', '%0.2g nM'));
legendStr2 = cellstr(num2str(mergedExpFile(8).concentrations{2}', '%0.2g nM'));

barweb(s70DataMtx,s70StdMtx,[],legendStr,'GFP protein expression endpoints ',[],[],'summer','y',legendStr2,[],'plot');

t_vec_mins=mergedExpFile(8).t_vec/60;

figure('Name','Pr-GFP +s70 kinetics');
hold on
h = stdshade(t_vec_mins,mergedExpFile(8).noBg_mean(:,:,2)/RFUumConvert,mergedExpFile(8).noBg_std(:,:,2)/RFUumConvert/sqrt(3),0.3,colorCodes);
axis([-10,850,-1000/RFUumConvert,35000/RFUumConvert])
xlabel('Time (mins)')
ylabel('deGFP fluorescence (A.U.)')
grid on
%legendStr = cellstr(num2str(mergedExpFile(8).concentrations(:)', 'Pr-deGFP-MGapt %0.2g nM'));
%legend(h,legendStr,'Location','NorthWest')

jj = reshape(1:36,3,12)';
p1 = reshape(jj(1:3:end,:)',1,12)+1;
p2 = reshape(jj(2:3:end,:)',1,12)+1;
p3 = reshape(jj(3:3:end,:)',1,12)+1;
%
% mergedExpFile(8) = mergeWholeExpFile({expFile(12)},[p1; p2; p3]+1);
% mergedExpFile(8).name = 'Pr-s70 pr-GFP-mg15';

%just use Pr-GFP as background for now:
mgBackground = repmat(mergedExpFile(1).Data_mean(:,[6 7 9],1),1,4);  %%  <- this is for 1,2,and 10 nM Pr-GFP

avgMGmean=mean([sum(expFile(12).Data(:,p1,1)-mgBackground); sum(expFile(12).Data(:,p2,1)-mgBackground); sum(expFile(12).Data(:,p3,1)-mgBackground)]);
avgMGerr=std([sum(expFile(12).Data(:,p1,1)-mgBackground); sum(expFile(12).Data(:,p2,1)-mgBackground); sum(expFile(12).Data(:,p3,1)-mgBackground)])/sqrt(3);

GFPendpoints=mergedExpFile(8).noBg_mean(end,1:12,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(8).noBg_std(end,1:12,2)/RFUumConvert/sqrt(3);

figure
clear legh
ylim([2 16])
xlim([50000 300000])

hold on
for k=0:3
    idx = [1:3]+k*3;
    legh(k+1) = hverrorbar(avgMGmean(idx),GFPendpoints(idx),avgMGerr(idx),GFPendpoints_err(idx),colorCodes{k+1});
end
%%%%%  TOTAL MG VS GFP
avgMGmean=mean([sum(expFile(3).Data(:,2:11,1)-mgNTPBackground); sum(expFile(3).Data(:,12:21,1)-mgNTPBackground); sum(expFile(3).Data(:,22:31,1)-mgNTPBackground)]);
avgMGerr=std([sum(expFile(3).Data(:,2:11,1)-mgNTPBackground); sum(expFile(3).Data(:,12:21,1)-mgNTPBackground); sum(expFile(3).Data(:,22:31,1)-mgNTPBackground)])/sqrt(3);

GFPendpoints=mergedExpFile(2).noBg_mean(end,:,2)/RFUumConvert;
GFPendpoints_err=mergedExpFile(2).noBg_std(end,:,2)/RFUumConvert/sqrt(3);

legh(end+1) = hverrorbar(avgMGmean,GFPendpoints,avgMGerr,GFPendpoints_err,'k');


legend(legh,'Pr-deGFP-MGApt + 0.1 nM Pr-\sigma70','Pr-deGFP-MGApt + 0.2 nM Pr-\sigma70',...
    'Pr-deGFP-MGApt + 0.5 nM Pr-\sigma70','Pr-deGFP-MGApt + 1 nM Pr-\sigma70',...
    'Pr-deGFP-MGApt','Location','SouthEast')

hold off
ylabel('deGFP (uM)')
xlabel('total MG r.f.u.')

grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr-gfp pr-CFP-mg15

GFP_DataMtx = reshape(mergedExpFile(7).noBg_mean(end,:,2),4,4)';
GFP_StdMtx  = reshape(mergedExpFile(7).noBg_std(end,:,2),4,4)';

CFP_DataMtx = reshape(mergedExpFile(7).noBg_mean(end,:,3),4,4)';
CFP_StdMtx  = reshape(mergedExpFile(7).noBg_std(end,:,3),4,4)';

% figure('Name','endpoint GFP')
% legendStr = cellstr(num2str(mergedExpFile(7).concentrations', '%0.2g nM'));
% barweb(GFP_DataMtx,GFP_StdMtx,[],legendStr,'GFP protein expression endpoints ',[],[],'summer','y',legendStr,[],'plot');
%
% figure('Name','endpoint CFP')
% barweb(CFP_DataMtx,CFP_StdMtx,[],legendStr,'CFP protein expression endpoints ',[],[],'summer','y',legendStr,[],'plot');

jj = reshape(1:32,4,8)';
p1 = jj(1:2,:);
p2 = jj(3:4,:);
p3 = jj(5:6,:);
p4 = jj(7:8,:);

q=[p1 p2 p3 p4; 1:16]+1;

mgBackground = repmat(mergedExpFile(1).Data_mean(:,[6 7 8 9],1),1,4);
totalTime=expFile(10).t_vec(end)/60;

avgMGmean=mean([sum(max(expFile(10).Data(:,q(1,:),1)-mgBackground,0)); sum(max(expFile(10).Data(:,q(2,:),1)-mgBackground,0)); sum(max(expFile(11).Data(:,q(3,:),1)-mgBackground,0))])*3/RFUumConvertMG/60;
avgMGerr=std([sum(max(expFile(10).Data(:,q(1,:),1)-mgBackground,0)); sum(max(expFile(10).Data(:,q(2,:),1)-mgBackground,0)); sum(max(expFile(11).Data(:,q(3,:),1)-mgBackground,0))])/sqrt(3)*3/RFUumConvertMG/60;

GFPendpoints    = mergedExpFile(7).noBg_mean(end,:,2)/RFUumConvert;
GFPendpointsstd = mergedExpFile(7).noBg_std(end,:,2)/RFUumConvert;

CFPendpoints    = mergedExpFile(7).noBg_mean(end,:,3);
CFPendpointsstd = mergedExpFile(7).noBg_std(end,:,3);

%%% begin figure %%%
figure('Name','cfp vs mg')
hold on
for k=0:3
    idx =  [1:4]+k*4;
    h(k+1) = hverrorbar(avgMGmean(idx),CFPendpoints(idx),avgMGerr(idx),CFPendpointsstd(idx),colorCodes{k+1});
end
hold off
xlabel('avg MG (nM)')
ylabel('deCFP (R.F.U.)')
legend(h,'pr-gfp 1-10nM + pr-cfp-MG 1nM','pr-gfp 1-10nM + pr-cfp-MG 2nM','pr-gfp 1-10nM + pr-cfp-MG 5nM','pr-gfp 1-10nM + pr-cfp-MG 10nM') ;
grid on
axis([0 120 0 14000])
%%% end figure %%%


%%% begin figure %%%
figure('Name','gfp vs cfp')
hold on
for k=0:3
    idx =  [1:4]+k*4;
    h(k+1) = hverrorbar(CFPendpoints(idx),GFPendpoints(idx),CFPendpointsstd(idx),GFPendpointsstd(idx),colorCodes{k+1});
end
hold off
xlabel('deCFP r.f.u.')
ylabel('deGFP (uM)')
legend(h,'pr-gfp 1-10nM + pr-cfp-MG 1nM','pr-gfp 1-10nM + pr-cfp-MG 2nM','pr-gfp 1-10nM + pr-cfp-MG 5nM','pr-gfp 1-10nM + pr-cfp-MG 10nM') ;
grid on

%%% end figure %%%


%%% begin figure %%%
figure('Name','gfp vs avg MG')
hold on

for k=0:3
    idx =  [1:4]+k*4;
    h(k+1) = hverrorbar(avgMGmean(idx),GFPendpoints(idx),avgMGerr(idx),GFPendpointsstd(idx),colorCodes{k+1});
end
hold off
xlabel('avg MG (nM)')
ylabel('deGFP (uM)')
legend(h,'pr-gfp 1-10nM + pr-cfp-MG 1nM','pr-gfp 1-10nM + pr-cfp-MG 2nM','pr-gfp 1-10nM + pr-cfp-MG 5nM','pr-gfp 1-10nM + pr-cfp-MG 10nM') ;
grid on
%%% end figure %%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pr-gfp_delRBS pr-cfp-mg15

GFP_DataMtx = reshape(mergedExpFile(9).noBg_mean(end,:,2),4,4)';
GFP_StdMtx  = reshape(mergedExpFile(9).noBg_std(end,:,2),4,4)';

CFP_DataMtx = reshape(mergedExpFile(9).noBg_mean(end,:,3),4,4)';
CFP_StdMtx  = reshape(mergedExpFile(9).noBg_std(end,:,3),4,4)';

% figure('Name','endpoint GFP')
% legendStr = cellstr(num2str(mergedExpFile(9).concentrations', '%0.2g nM'));
% barweb(GFP_DataMtx,GFP_StdMtx,[],legendStr,'GFP protein expression endpoints ',[],[],'summer','y',legendStr,[],'plot');
%
% figure('Name','endpoint CFP')
% barweb(CFP_DataMtx,CFP_StdMtx,[],legendStr,'CFP protein expression endpoints ',[],[],'summer','y',legendStr,[],'plot');

jj = reshape(1:48,4,12)';
p1 = jj(1:3,:);
p2 = jj(4:6,:);
p3 = jj(7:9,:);
p4 = jj(10:12,:);

q=[p1 p2 p3 p4]+1;

mgBackground = repmat(mergedExpFile(1).Data_mean(:,[6 7 8 9],1),1,4);
totalTime=expFile(13).t_vec(end)/60;

avgMGmean=mean([sum(max(expFile(13).Data(:,q(1,:),1)-mgBackground,0)); sum(max(expFile(13).Data(:,q(2,:),1)-mgBackground,0)); sum(max(expFile(13).Data(:,q(3,:),1)-mgBackground,0))])*3/RFUumConvertMG/60;
avgMGerr=std([sum(max(expFile(13).Data(:,q(1,:),1)-mgBackground,0)); sum(max(expFile(13).Data(:,q(2,:),1)-mgBackground,0)); sum(max(expFile(13).Data(:,q(3,:),1)-mgBackground,0))])/sqrt(3)*3/RFUumConvertMG/60;

GFPendpoints    = mergedExpFile(9).noBg_mean(end,:,2)/RFUumConvert;
GFPendpointsstd = mergedExpFile(9).noBg_std(end,:,2)/RFUumConvert;

CFPendpoints    = mergedExpFile(9).noBg_mean(end,:,3);
CFPendpointsstd = mergedExpFile(9).noBg_std(end,:,3);

%%% begin figure %%%
figure('Name','cfp vs mg')
hold on
for k=0:3
    idx =  [1:4]+k*4;
    h(k+1) = hverrorbar(avgMGmean(idx),CFPendpoints(idx),avgMGerr(idx),CFPendpointsstd(idx),colorCodes{k+1});
end
hold off
xlabel('avg MG (nM)')
ylabel('deCFP (R.F.U.)')
legend(h,'pr-gfp 1-10nM + pr-cfp-MG 1nM','pr-gfp 1-10nM + pr-cfp-MG 2nM','pr-gfp 1-10nM + pr-cfp-MG 5nM','pr-gfp 1-10nM + pr-cfp-MG 10nM') ;
grid on
axis([0 120 0 14000])
%%% end figure %%%


%%% begin figure %%%
figure('Name','gfp vs cfp')
hold on
for k=0:3
    idx =  [1:4]+k*4;
    h(k+1) = hverrorbar(CFPendpoints(idx),GFPendpoints(idx),CFPendpointsstd(idx),GFPendpointsstd(idx),colorCodes{k+1});
end
hold off
xlabel('deCFP r.f.u.')
ylabel('deGFP (uM)')
legend(h,'pr-gfp 1-10nM + pr-cfp-MG 1nM','pr-gfp 1-10nM + pr-cfp-MG 2nM','pr-gfp 1-10nM + pr-cfp-MG 5nM','pr-gfp 1-10nM + pr-cfp-MG 10nM') ;
grid on
%%% end figure %%%


%%% begin figure %%%
figure('Name','gfp vs avg MG')
hold on

for k=0:3
    idx =  [1:4]+k*4;
    h(k+1) = hverrorbar(avgMGmean(idx),GFPendpoints(idx),avgMGerr(idx),GFPendpointsstd(idx),colorCodes{k+1});
end
hold off
xlabel('avg MG (nM)')
ylabel('deGFP (uM)')
legend(h,'pr-gfp 1-10nM + pr-cfp-MG 1nM','pr-gfp 1-10nM + pr-cfp-MG 2nM','pr-gfp 1-10nM + pr-cfp-MG 5nM','pr-gfp 1-10nM + pr-cfp-MG 10nM') ;
grid on
%%% end figure %%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pr-Cfp-mg15


dataMtx = reshape(mergedExpFile(10).noBg_mean(end,:,2),1,10);
stdMtx  = reshape(mergedExpFile(10).noBg_std(end,:,2),1,10);


small_amounts = 1:5;
large_amounts = 6:10;
figure('Name','endpoint small')
barweb(dataMtx(small_amounts),stdMtx(small_amounts),[],mergedExpFile(10).concentrations(small_amounts),'CFP protein expression endpoints ',[],[],'summer','y',mergedExpFile(10).consturctNames,[],'plot');
figure('Name','endpoint large')
barweb(dataMtx(large_amounts),stdMtx(large_amounts),[],mergedExpFile(10).concentrations(large_amounts),'CFP protein expression endpoints ',[],[],'summer','y',mergedExpFile(10).consturctNames,[],'plot');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


