data_loc = 'C:\Users\ZappaLab\Google Drive\Experiments\2020 SST Rain Response Paper\Data\Working\   ';

load([data_loc '\SST_rain_comparison_peak_data.mat'])
load([data_loc '\SST_rain_comparison_variables2.mat'])

% ----------------------------------------------------
%  Some basic plots
% ----------------------------------------------------

%may need to change these to column vectors at some point
SST_max = [Vars.Max_Response]; 
SST_Tmax = [Vars.Time_to_Max];
Rain_Length = [Vars.Length_of_Rain_Event];
Cum_Rain_Rate = [Vars.Total_Rain];
SST0 = [Vars.SST0];
Max_Response = SST0-SST_max;
Max_Rain = [Vars.Peak];
Time_to_Peak_Rain = Max_Rain*NaN;
npeaks = size(Vars,2);


for i = 1:npeaks
   Time_to_Peak_Rain(i) = Vars(i).Time_of_Peak - Vars(i).Time(1+beg_pad1);
   Time_to_Peak_Rain(i) = Time_to_Peak_Rain(i) * 24 * 60;
end



figure(1);clf
orient landscape

subplot(2,3,1)
bins = 0:0.05:2.6;
N = histc((Max_Response),bins );
%plot(bins(1:end-1)+0.5*(diff(bins(1:2))), N(1:end-1),'.-')
bar(bins,N,'histc')
xlim([0 2.7])
set(gca,'FontSize',12)
title('Max \DeltaSST_{skin}','FontSize',14)
xlabel('[\circC]','FontSize',14)

subplot(2,3,2)
bins = 0:0.05:1;
N = histc((SST_Tmax./Rain_Length),bins );
%plot(bins(1:end-1)+0.5*(diff(bins(1:2))), N(1:end-1),'.-')
bar(bins,N,'histc')
title('t_{*SST}','FontSize',14)
xlabel('((t-t_0)/L_{rain}','FontSize',14)
%xlabel('[Minutes]','FontSize',14)
ylim([0 15])
set(gca,'FontSize',12)
xlim([0 1.1])


subplot(2,3,3)
bins = (0:10:300); 
N = histc(SST_Tmax,bins);
%plot(bins(1:end-1)+0.5*(diff(bins(1:2))), N(1:end-1),'.-')
bar(bins,N,'histc')
title('t_{*SST}','FontSize',14)
xlabel('[minutes]','FontSize',14)
set(gca,'FontSize',12)
xlim([0 250])

subplot(2,3,4)
bins = 0:1:900;
N = histc(Cum_Rain_Rate,bins);
%plot(bins(1:end-1)+0.5*(diff(bins(1:2))), N(1:end-1),'.-')
bar(bins,N,'histc')
title('Cumulative Rain','FontSize',14)
xlabel('[mm]','FontSize',14)
xlim([0 150])
set(gca,'FontSize',12)

subplot(2,3,5)
bins = 0:0.25:15;
N = histc(Rain_Length/60,bins);
%plot(bins(1:end-1)+0.5*(diff(bins(1:2))), N(1:end-1),'.-')
bar(bins,N,'histc')
title('L_{rain}','FontSize',14)
xlabel('[hours]','FontSize',14)
set(gca,'FontSize',12)
xlim([0 16])


subplot(2,3,6)
bins = 0:15:475;
N = histc(Time_to_Peak_Rain,bins);
bar(bins,N,'histc')
title('t_{*rain}','FontSize',14)
xlabel('[minutes]','FontSize',14)
xlim([0 250])
set(gca,'FontSize',12)


%print -dpng Dynamo_RainEvent_histograms.png
%print -dpsc Dynamo_RainEvent_histograms.ps

%orient portrait
%set(gcf,'papersize',[11 8.5])
%set(gcf,'paperposition',[.25 .25 10.5 8])
%print -dpng Dynamo_RainEvent_histograms.png

%return





figure(2);clf
orient landscape

subplot(2,3,1)
plot((SST0-SST_max), Cum_Rain_Rate,'.')
xlabel('Max SST_{skin} Response [\circC]','FontSize',14)
ylabel('Cumulative Rain [mm]','FontSize',14)
ylim([0 400])

subplot(2,3,4)
plot((SST_Tmax./Rain_Length), Cum_Rain_Rate,'.')
ylabel('Cumulative Rain [mm]','FontSize',14)
xlabel('Time to max response/event duration','FontSize',14)
ylim([0 400])


Max_Rain = SST0*NaN;
WS = SST0*NaN;
for i = 1:size(Vars,2)
   ind = 1:length(Vars(i).Rain);
   ind(1:2) = []; 
   ind(end-3:end) = [];
   tmp = Vars(i).Rain;
   Max_Rain(i) = max(tmp(ind)); 
   WS(i) = Vars(i).Wind(3);
end



subplot(2,3,2)
plot((SST0-SST_max), Max_Rain,'.')
xlabel('Max SST_{skin} Response [\circC]','FontSize',14)
ylabel('Max Rain Rate [mm]','FontSize',14)

subplot(2,3,5)
plot((SST_Tmax./Rain_Length), Max_Rain,'.')
ylabel('Max Rain Rate [mm]','FontSize',14)
xlabel('Time to max response/event duration','FontSize',14)

subplot(2,3,3)
plot((SST0-SST_max), Rain_Length,'.')
xlabel('Max SST_{skin} Response [\circC]','FontSize',14)
ylabel('Length of Rain Event [minutes]','FontSize',14)
ylim([0 600])

subplot(2,3,6)
plot((SST_Tmax./Rain_Length), Rain_Length,'.')
ylabel('Length of Rain Event [minutes]','FontSize',14)
xlabel('Time to max response/event duration','FontSize',14)
ylim([0 600])

WSlow = find(WS < 5);
WShigh = find(WS > 9);
WSmed = find(WS <= 9 | WS >= 5);

%print -dpng Dynamo_Rain_scatterplots.png
%print -dpsc Dynamo_Rain_scatterplots.ps




figure(3);clf
orient landscape

subplot(2,4,1)
plot((SST0-SST_max), Cum_Rain_Rate,'.')
xlabel('Max \DeltaSST_{skin} [\circC]','FontSize',14)
ylabel('Cumulative Rain [mm]','FontSize',14)
ylim([0 400])
set(gca,'FontSize',12)

subplot(2,4,5)
plot((SST_Tmax), Cum_Rain_Rate,'.')
ylabel('Cumulative Rain [mm]','FontSize',14)
xlabel('t_{*SST} [minutes]','FontSize',14)
ylim([0 400])
set(gca,'FontSize',12)

%add_datestamp


Max_Rain = SST0*NaN;
WS = SST0*NaN;
for i = 1:size(Vars,2)
   ind = 1:length(Vars(i).Rain);
   ind(1:2) = [];
   ind(end-3:end) = [];
   tmp = Vars(i).Rain;
   Max_Rain(i) = max(tmp(ind));
   WS(i) = Vars(i).Wind(3);
end



subplot(2,4,2)
plot((SST0-SST_max), Max_Rain,'.')
xlabel('Max \DeltaSST_{skin} [\circC]','FontSize',14)
ylabel('Max Rain Rate [mm]','FontSize',14)
set(gca,'FontSize',12)

subplot(2,4,6)
plot((SST_Tmax), Max_Rain,'.')
ylabel('Max Rain Rate [mm]','FontSize',14)
xlabel('t_{*SST} [minutes]','FontSize',14)
set(gca,'FontSize',12)

subplot(2,4,3)
plot((SST0-SST_max), Rain_Length,'.')
xlabel('Max \DeltaSST_{skin} [\circC]','FontSize',14)
ylabel('L_{rain} [minutes]','FontSize',14)
ylim([0 600])
set(gca,'FontSize',12)


subplot(2,4,7)
plot((SST_Tmax), Rain_Length,'.')
ylabel('L_{rain} [minutes]','FontSize',14)
xlabel('t_{*SST} [minutes]','FontSize',14)
ylim([0 600])
hold on
plot([0 400],[0 400],'linewidth',1)
set(gca,'FontSize',12)



subplot(2,4,4)
plot((SST0-SST_max), Rain_Length,'.')
xlabel('Max \DeltaSST_{skin} [\circC]','FontSize',14)
ylabel('t_{*rain} [minutes]','FontSize',14)
%ylim([0 600])
set(gca,'FontSize',12)

subplot(2,4,8)
plot(SST_Tmax,Time_to_Peak_Rain,'.')
ylabel('t_{*rain} [minutes]','FontSize',14)
xlabel('t_{*SST} [minutes]','FontSize',14)
hold on
plot([0 500],[0 500],'linewidth',1)
set(gca,'FontSize',12)




%print -dpng Dynamo_Rain_scatterplots_unscaled.png
%print -dpsc Dynamo_Rain_scatterplots_unscaled.ps
%saveas(gcf,'Dynamo_Rain_scatterplots_unscaled.fig','fig')
%orient portrait
%set(gcf,'papersize',[11 8.5])
%set(gcf,'paperposition',[.25 .25 10.5 8])
%print -dpng Dynamo_Rain_scatterplots_unscaled.png



%figure(4);clf
%orient landscape
%subplot(2,2,1)
%plot([Vars.SST0]-[Vars.Max_Response],-1*[Vars.Tair],'o','MarkerSize',4,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
%ylabel('\DeltaSST [\circC]','FontSize',14)
%ylabel('\DeltaT_{air}[\circC]','FontSize',14)

%RHFsum = zeros(npeaks,1);
%for i = 1:npeaks
%   RHFsum(i) = sum(Vars(i).RHF);
%end

%subplot(2,2,2)
%plot(RHFsum,col2(Vars.Tair_response),'o','MarkerSize',4,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
%semilogx(-1*RHFsum,[Vars.Tair_response],'o','MarkerSize',4,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
%ylabel('\DeltaT_{air}[\circC]','FontSize',14)
%xlabel('Latent Heat Flux Due to Rain','FontSize',14)

%subplot(2,2,3)
%semilogx(RHFsum,[Vars.SST0]-[Vars.Max_Response],'o','MarkerSize',4,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
%ylabel('\DeltaSST [\circC]','FontSize',14)
%xlabel('Latent Heat Flux Due to Rain','FontSize',14)


%print -dpng Dynamo_Rain_scatterplots_unscaled2.png
%print -dpsc Dynamo_Rain_scatterplots_unscaled2.ps
%saveas(gcf,'Dynamo_Rain_scatterplots_unscaled2.fig','fig')





%return

%WSlow = find(WS < 5);
%WShigh = find(WS > 9);
%WSmed = find(WS <= 9 | WS >= 5);









figure(5);clf
orient landscape

subplot(2,3,1)
plot((SST0(WSlow)-SST_max(WSlow)), Cum_Rain_Rate(WSlow),'.')
hold on
plot((SST0(WSmed)-SST_max(WSmed)), Cum_Rain_Rate(WSmed),'g.')
plot((SST0(WShigh)-SST_max(WShigh)), Cum_Rain_Rate(WShigh),'r.')
xlabel('Max SST_{skin} Response [\circC]','FontSize',14)
ylabel('Cumulative Rain [mm]','FontSize',14)

subplot(2,3,4)
plot((SST_Tmax(WSlow)./Rain_Length(WSlow)), Cum_Rain_Rate(WSlow),'.')
hold on
plot((SST_Tmax(WSmed)./Rain_Length(WSmed)), Cum_Rain_Rate(WSmed),'g.')
plot((SST_Tmax(WShigh)./Rain_Length(WShigh)), Cum_Rain_Rate(WShigh),'r.')
ylabel('Cumulative Rain [mm]','FontSize',14)
xlabel('Time to max response/event duration','FontSize',14)

subplot(2,3,2)
plot((SST0(WSlow)-SST_max(WSlow)), Max_Rain(WSlow),'.')
hold on
plot((SST0(WSmed)-SST_max(WSmed)), Max_Rain(WSmed),'g.')
plot((SST0(WShigh)-SST_max(WShigh)), Max_Rain(WShigh),'r.')
xlabel('Max SST_{skin} Response [\circC]','FontSize',14)
ylabel('Max Rain Rate [mm]','FontSize',14)

subplot(2,3,5)
plot((SST_Tmax(WSlow)./Rain_Length(WSlow)), Max_Rain(WSlow),'.')
hold on
plot((SST_Tmax(WSmed)./Rain_Length(WSmed)), Max_Rain(WSmed),'g.')
plot((SST_Tmax(WShigh)./Rain_Length(WShigh)), Max_Rain(WShigh),'r.')
ylabel('Max Rain Rate [mm]','FontSize',14)
xlabel('Time to max response/event duration','FontSize',14)

subplot(2,3,3)
plot((SST0(WSlow)-SST_max(WSlow)), Rain_Length(WSlow),'.')
hold on
plot((SST0(WSmed)-SST_max(WSmed)), Rain_Length(WSmed),'g.')
plot((SST0(WShigh)-SST_max(WShigh)), Rain_Length(WShigh),'r.')
xlabel('Max SST_{skin} Response [\circC]','FontSize',14)
ylabel('Length of Rain Event [minutes]','FontSize',14)

subplot(2,3,6)
plot((SST_Tmax(WSlow)./Rain_Length(WSlow)), Rain_Length(WSlow),'.')
hold on
plot((SST_Tmax(WSmed)./Rain_Length(WSmed)), Rain_Length(WSmed),'g.')
plot((SST_Tmax(WShigh)./Rain_Length(WShigh)), Rain_Length(WShigh),'r.')

ylabel('Length of Rain Event [minutes]','FontSize',14)
xlabel('Time to max response/event duration','FontSize',14)








%-----------------------------------------------------
%   Bin average rain events
%-----------------------------------------------------

maxL = -32768;
minL = 32768;
npeaks = size(Vars,2);
dt1 = 1/24;		% one hour gap in Edson's data  

for i = 1:npeaks
   if length(Vars(i).TRel21) > 0
      if Vars(i).TRel21(end)*24*60./Vars(i).Time_to_Max > maxL
         maxL = Vars(i).TRel21(end)*24*60./Vars(i).Time_to_Max;
      end
      if Vars(i).TRel21(1)*24*60./Vars(i).Time_to_Max < minL
         minL = Vars(i).TRel21(1)*24*60./Vars(i).Time_to_Max;
      end
   end
end


C = jet(npeaks);



% Do some bin-averaging

%bins = 0.9989:0.000025:1.00120;
bins = floor(minL):0.5:ceil(maxL);
bins = bins + 0.5*(diff(bins(1:2)));

nmnSST = zeros(length(bins),1);
nsdSST = zeros(length(bins),1);
nmnbulk = zeros(length(bins),1);
nsdbulk = zeros(length(bins),1);
nmnTair = zeros(length(bins),1);
nsdTair = zeros(length(bins),1);
nmnTbulk = zeros(length(bins),1);
nsdTbulk = zeros(length(bins),1);
nmnP = zeros(length(bins),1);
nsdP = zeros(length(bins),1);
nmnU10r = zeros(length(bins),1);
nsdU10r = zeros(length(bins),1);
npnts = zeros(length(bins),1);
npnts2 = zeros(length(bins),1);


% Try using change from initial condition rather than straight value
% calculate mean values in bins

for i = 1:length(bins)
   Bins(i).Tair = NaN;
   Bins(i).U10r = NaN;
   Bins(i).Tbulk = NaN;
   Bins(i).SST = NaN;
   Bins(i).T1 = NaN;
end

for j = 1:npeaks
   %if ~isempty(Vars(j).Time_to_Max)
   if ~isnan(Vars(j).Time_to_Max)
      TRel1 = Vars(j).TRel11*24*60./Vars(j).Time_to_Max;
      Tair = Vars(j).Tair;
      U10r = Vars(j).Wind;
      Tbulk = Vars(j).Tbulk;
      bad = find(diff(Vars(j).Time) > dt1);
      TRel1(bad+1:end) = [];
      Tair(bad+1:end) = [];
      U10r(bad+1:end) = [];
      Tbulk(bad+1:end) = [];

      Tair = Tair - Tair(1+beg_pad1);
      U10r = U10r - U10r(1+beg_pad1);
      Tbulk = Tbulk - Tbulk(1+beg_pad1);
   
      nn = Vars(j).time_indices2;
      TRel2 = Vars(j).TRel21*24*60./Vars(j).Time_to_Max;
      SST = Vars(j).SST_smoothed;
      indx1 = 1+beg_pad2:-1:1+beg_pad2-18;
      SST = SST - nanmean(SST(indx1));

      indx2 = NaN*ones(length(Tair),1);
      for k = 1:length(TRel1)
         indx = find(bins > TRel1(k) );
         indx = indx(1)-1;
         tmp = Bins(indx).U10r;
         tmp = [tmp U10r(k)];
         Bins(indx).U10r = tmp;
         tmp = Bins(indx).Tair;
         tmp = [tmp Tair(k)];
         Bins(indx).Tair= tmp;
         tmp = Bins(indx).Tbulk;
         tmp = [tmp Tbulk(k)];
         Bins(indx).Tbulk = tmp;
         tmp = Bins(indx).T1;
         tmp = [tmp TRel1(k)];
         Bins(indx).T1 = tmp;
      end
      for k = 1:length(TRel2)
         indx = find(bins > TRel2(k) );
         indx = indx(1)-1;
         %nmnSST(indx) = nmnSST(indx) + SST(k);
         %npnts2(indx) = npnts2(indx) + 1;
         tmp = Bins(indx).SST;
         tmp = [tmp SST(k)];
         Bins(indx).SST = tmp;
      end

   end
end


for i = 1:length(bins)
   nmnTair(i) = nanmean(Bins(i).Tair);
   nsdTair(i) = nanstd(Bins(i).Tair);
   nmnTbulk(i) = nanmean(Bins(i).Tbulk);
   nsdTbulk(i) = nanstd(Bins(i).Tbulk);
   nmnU10r(i) = nanmean(Bins(i).U10r);
   nsdU10r(i) = nanstd(Bins(i).U10r);
   nmnSST(i) = nanmean(Bins(i).SST);
   nsdSST(i) = nanstd(Bins(i).SST);
end


figure(6);clf
orient tall
ax(1) = subplot(4,1,1);
plot(bins,nmnTair,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTair(i)-nsdTair(i) nmnTair(i)+nsdTair(i)],'k');
end
title('Dynamo Rain Events','FontSize',16)
ylabel('\DeltaT_{air} [\circC]','FontSize',14)
set(gca,'FontSize',12)
ylim([-3 1])
line([-3, 10],[0,0])


ax(2) = subplot(4,1,2);
plot(bins,nmnU10r,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnU10r(i)-nsdU10r(i) nmnU10r(i)+nsdU10r(i)],'k');
end
ylabel('\DeltaU_{r_{10}} [m/s]','FontSize',14)
set(gca,'FontSize',12)
ylim([-4 4])
line([-3, 10],[0,0])

ax(3) = subplot(4,1,3);
plot(bins,nmnTbulk,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTbulk(i)-nsdTbulk(i) nmnTbulk(i)+nsdTbulk(i)],'k');
end
ylabel('\DeltaT_{bulk} [\circC]','FontSize',14)
set(gca,'FontSize',12)
ylim([-0.5 0.2])
line([-3, 10],[0,0])


ax(4) = subplot(4,1,4);
plot(bins,nmnSST,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnSST(i)-nsdSST(i) nmnSST(i)+nsdSST(i)],'k');
end
ylabel('\DeltaSST_{skin} [\circC]','FontSize',14)
xlabel('(T-T_0)/t_{*SST}','FontSize',14)
ylim([-0.8 0.4])
line([-3, 10],[0,0])
set(gca,'FontSize',12,'ytick',[-0.8:0.4:0.4])
linkaxes(ax,'x')
xlim([-15 40])

%print -dpng Dynamo_RainEvent_Composite.png
%print -dpsc Dynamo_RainEvent_Composite.ps

set(ax(2),'xlim',[-3 10])

%saveas(gcf,'Dynamo_RainEvent_Composite_zoom.fig','fig')
%print -dpng Dynamo_RainEvent_Composite_zoom.png
%print -dpsc Dynamo_RainEvent_Composite_zoom.ps








%-----------------------------------------------------
%   Bin average rain events - do not scale time
%-----------------------------------------------------

maxL = -32768;
minL = 32768;
npeaks = size(Vars,2);
dt1 = 1/24;		% one hour gap in Edson's data  

for i = 1:npeaks
   if length(Vars(i).TRel21) > 0
      if Vars(i).TRel21(end)*24*60 > maxL
         maxL = Vars(i).TRel21(end)*24*60.;
      end
      if Vars(i).TRel21(1)*24*60 < minL
         minL = Vars(i).TRel21(1)*24*60;
      end
   end
end


C = jet(npeaks);



% Do some bin-averaging

%bins = 0.9989:0.000025:1.00120;
bins = floor(minL):10:ceil(maxL);

nmnSST = zeros(length(bins),1);
nsdSST = zeros(length(bins),1);
nmnbulk = zeros(length(bins),1);
nsdbulk = zeros(length(bins),1);
nmnTair = zeros(length(bins),1);
nsdTair = zeros(length(bins),1);
nmnTbulk = zeros(length(bins),1);
nsdTbulk = zeros(length(bins),1);
nmnP = zeros(length(bins),1);
nsdP = zeros(length(bins),1);
nmnU10r = zeros(length(bins),1);
nsdU10r = zeros(length(bins),1);
npnts = zeros(length(bins),1);
npnts2 = zeros(length(bins),1);


% Try using change from initial condition rather than straight value
% calculate mean values in bins

for j = 1:npeaks
   Bins(i).Tair = NaN;
   Bins(i).U10r = NaN;
   Bins(i).Tbulk = NaN;
   Bins(i).SST = NaN;
end

for j = 1:npeaks
   %if ~isempty(Vars(j).Time_to_Max)
   if ~isnan(Vars(j).Time_to_Max)
      TRel1 = Vars(j).TRel11*24*60;
      Tair = Vars(j).Tair;
      U10r = Vars(j).Wind;
      Tbulk = Vars(j).Tbulk;
      bad = find(diff(Vars(j).Time) > dt1);
      TRel1(bad+1:end) = [];
      Tair(bad+1:end) = [];
      U10r(bad+1:end) = [];
      Tbulk(bad+1:end) = [];

      Tair = Tair - Tair(1+beg_pad1);
      U10r = U10r - U10r(1+beg_pad1);
      Tbulk = Tbulk - Tbulk(1+beg_pad1);
   
      nn = Vars(j).time_indices2;
      TRel2 = Vars(j).TRel21*24*60;
      SST = Vars(j).SST_smoothed;
      indx1 = 1+beg_pad2:-1:1+beg_pad2-18;
      SST = SST - nanmean(SST(indx1));

      indx2 = NaN*ones(length(Tair),1);
      for k = 1:length(TRel1)
         indx = find(bins > TRel1(k) );
         indx = indx(1)-1;
         tmp = Bins(indx).U10r;
         tmp = [tmp U10r(k)];
         Bins(indx).U10r = tmp;
         tmp = Bins(indx).Tair;
         tmp = [tmp Tair(k)];
         Bins(indx).Tair= tmp;
         tmp = Bins(indx).Tbulk;
         tmp = [tmp Tbulk(k)];
         Bins(indx).Tbulk = tmp;

         tmp = Bins(indx).SST;
         tmp = [tmp SST(k)];
         Bins(indx).SST= tmp;

      end
      for k = 1:length(TRel2)
         indx = find(bins > TRel2(k) );
         indx = indx(1)-1;
         nmnSST(indx) = nmnSST(indx) + SST(k); 
         npnts2(indx) = npnts2(indx) + 1;
         tmp = Bins(indx).SST;
         tmp = [tmp SST(k)];
         Bins(indx).SST= tmp;
      end
   end
end

for i = 1:length(bins)
   nmnTair(i) = nanmean(Bins(i).Tair);
   nsdTair(i) = nanstd(Bins(i).Tair);
   nmnTbulk(i) = nanmean(Bins(i).Tbulk);
   nsdTbulk(i) = nanstd(Bins(i).Tbulk);
   nmnU10r(i) = nanmean(Bins(i).U10r);
   nsdU10r(i) = nanstd(Bins(i).U10r);
   nmnSST(i) = nanmean(Bins(i).SST);
   nsdSST(i) = nanstd(Bins(i).SST);
end





figure(7);clf
orient tall
ax(1) = subplot(4,1,1)
plot(bins,nmnTair,'-o','MarkerSize',5')
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTair(i)-nsdTair(i) nmnTair(i)+nsdTair(i)]);
end
title('Dynamo Rain Events','FontSize',14)
ylabel('Tair [\circC]','FontSize',14)
ylim([-3 2])
%yline(0)


ax(2) = subplot(4,1,2)
plot(bins,nmnU10r,'-o','MarkerSize',5')
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnU10r(i)-nsdU10r(i) nmnU10r(i)+nsdU10r(i)]);
end
ylabel('U_r_{10} [m/s]','FontSize',14)
ylim([-6 6])
%yline(0)

ax(3) = subplot(4,1,3)
plot(bins,nmnTbulk,'-o','MarkerSize',5')
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTbulk(i)-nsdTbulk(i) nmnTbulk(i)+nsdTbulk(i)]);
end
ylabel('T_{bulk} [\circC]','FontSize',14)
%ylim([-0.6 0.4])
%yline(0)


ax(4) = subplot(4,1,4)
plot(bins,nmnSST,'-o','MarkerSize',5')
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnSST(i)-nsdSST(i) nmnSST(i)+nsdSST(i)]);
end
ylabel('SST_{skin} [\circC]','FontSize',14)
xlabel('Time from Event Initiation [minutes]','FontSize',14)
ylim([-0.8 0.4])
%yline(0)
linkaxes(ax,'x')
xlim([-50 950])

%print -dpng Dynamo_RainEvent_Composite_time_unscaled.png
%print -dpsc Dynamo_RainEvent_Composite_time_unscaled.ps

%set(ax(2),'xlim',[-50 310])

%print -dpng Dynamo_RainEvent_Composite_time_unscaled_zoom.png
%print -dpsc Dynamo_RainEvent_Composite_time_unscaled_zoom.ps









