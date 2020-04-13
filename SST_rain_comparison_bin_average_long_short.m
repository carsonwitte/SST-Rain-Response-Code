% need to hand-edit final following and first leading minima for each segment
% peakdet.m does not pick them up for some reason


load SST_rain_comparison_peak_data.mat
load  SST_rain_comparison_variables2.mat

loc = '/local/data/deadshot1/Analysis/DYNAMO/KT15/IRSST/';

anomaly = [1 6 14 32 33 48 51 64 87];
long_events = [40 74 82 99];				% > 5 hours
high_rain_rate = [29 31 39 40 42 59 66 72 78 85 92];	% > 35 mm/hr
high_wind = [47 74 75 76 78 80 83 85 99 100];		% > 12 m/s
weak_events = [8 9 20 32 41 48 51 56 60 73 87];         % < 0.8 mm/hr




%-----------------------------------------------------
%   Bin average rain events  - Short only
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
      if ismember(j,long_events) | ismember(j,anomaly)
         continue
      end
      TRel1 = Vars(j).TRel11*24*60./Vars(j).Time_to_Max;
      Tair = Vars(j).Tair;
      U10r = Vars(j).Wind;
      Tbulk = Vars(j).Tbulk;
      Rain = Vars(j).Rain;
      Tme = Vars(j).Time;

      Tair = Tair - Tair(1+beg_pad1);
      U10r = U10r - U10r(1+beg_pad1);
      Tbulk = Tbulk - Tbulk(1+beg_pad1);
   
      nn = Vars(j).time_indices2;
      TRel2 = Vars(j).TRel21*24*60./Vars(j).Time_to_Max;
      SST = Vars(j).SST_smoothed;
      indx1 = 1+beg_pad2:-1:1+beg_pad2-18;
      SST = SST - nmean(SST(indx1));

      % Remove subsequent rain event, if present
      oo = Rain(end-beg_pad1:end);
      ooind = length(Rain)-beg_pad1:length(Rain);
      bad = find(oo > 0);
      if ~isempty(bad)
         Tbulk(ooind(bad)) = [];
         Tair(ooind(bad)) = [];
         TRel1(ooind(bad)) = [];
         U10r(ooind(bad)) = [];
         Tme(ooind(bad)) = [];
         bad2 = find(Vars(j).T2 > Tme(end));
         SST(bad2) = [];
         TRel2(bad2) = [];
      end
      % Remove data after large gaps
      bad = find(diff(Vars(j).Time) > dt1);
      TRel1(bad+1:end) = [];
      Tair(bad+1:end) = [];
      U10r(bad+1:end) = [];
      Tbulk(bad+1:end) = [];
      Rain(bad+1:end) = [];
      Tme(bad+1:end) = [];


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
   nmnTair(i) = nmean(Bins(i).Tair);
   nsdTair(i) = nstd(Bins(i).Tair);
   nmnTbulk(i) = nmean(Bins(i).Tbulk);
   nsdTbulk(i) = nstd(Bins(i).Tbulk);
   nmnU10r(i) = nmean(Bins(i).U10r);
   nsdU10r(i) = nstd(Bins(i).U10r);
   nmnSST(i) = nmean(Bins(i).SST);
   nsdSST(i) = nstd(Bins(i).SST);
end







figure
orient tall
ax(1) = subplot(4,1,1)
plot(bins,nmnTair,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTair(i)-nsdTair(i) nmnTair(i)+nsdTair(i)],'k');
end
title('Dynamo Rain Events - Events > 5 Hours Excluded','FontSize',16)
ylabel('\DeltaT_{air} [\circC]','FontSize',14)
set(gca,'FontSize',12)
ylim([-3 1])
hline(0)


ax(2) = subplot(4,1,2)
plot(bins,nmnU10r,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnU10r(i)-nsdU10r(i) nmnU10r(i)+nsdU10r(i)],'k');
end
ylabel('\DeltaU_{r_{10}} [m/s]','FontSize',14)
set(gca,'FontSize',12)
ylim([-4 4])
hline(0)

ax(3) = subplot(4,1,3)
plot(bins,nmnTbulk,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTbulk(i)-nsdTbulk(i) nmnTbulk(i)+nsdTbulk(i)],'k');
end
ylabel('\DeltaT_{bulk} [\circC]','FontSize',14)
set(gca,'FontSize',12)
ylim([-0.5 0.2])
hline(0)


ax(4) = subplot(4,1,4)
plot(bins,nmnSST,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnSST(i)-nsdSST(i) nmnSST(i)+nsdSST(i)],'k');
end
ylabel('\DeltaSST_{skin} [\circC]','FontSize',14)
xlabel('(T-T_0)/t_{*SST}','FontSize',14)
ylim([-0.8 0.4])
hline(0)
set(gca,'FontSize',12,'ytick',[-0.8:0.4:0.4])
linkaxes(ax,'x')
xlim([-15 40])

print -dpng Dynamo_RainEvent_Composite_short.png
print -dpsc Dynamo_RainEvent_Composite_short.ps

set(ax(2),'xlim',[-3 10])

saveas(gcf,'Dynamo_RainEvent_Composite_short_zoom.fig','fig')
print -dpng Dynamo_RainEvent_Composite_short_zoom.png
print -dpsc Dynamo_RainEvent_Composite_short_zoom.ps


% ---------------------
% Long Events
% ---------------------




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
      if ~ismember(j,long_events)
         continue
      end
      TRel1 = Vars(j).TRel11*24*60./Vars(j).Time_to_Max;
      Tair = Vars(j).Tair;
      U10r = Vars(j).Wind;
      Tbulk = Vars(j).Tbulk;
      Rain = Vars(j).Rain;
      Tme = Vars(j).Time;

      Tair = Tair - Tair(1+beg_pad1);
      U10r = U10r - U10r(1+beg_pad1);
      Tbulk = Tbulk - Tbulk(1+beg_pad1);
   
      nn = Vars(j).time_indices2;
      TRel2 = Vars(j).TRel21*24*60./Vars(j).Time_to_Max;
      SST = Vars(j).SST_smoothed;
      indx1 = 1+beg_pad2:-1:1+beg_pad2-18;
      SST = SST - nmean(SST(indx1));

      % Remove subsequent rain event, if present
      oo = Rain(end-beg_pad1:end);
      ooind = length(Rain)-beg_pad1:length(Rain);
      bad = find(oo > 0);
      if ~isempty(bad)
         Tbulk(ooind(bad)) = [];
         Tair(ooind(bad)) = [];
         TRel1(ooind(bad)) = [];
         U10r(ooind(bad)) = [];
         Tme(ooind(bad)) = [];
         bad2 = find(Vars(j).T2 > Tme(end));
         SST(bad2) = [];
         TRel2(bad2) = [];
      end
      % Remove data after large gaps
      bad = find(diff(Vars(j).Time) > dt1);
      TRel1(bad+1:end) = [];
      Tair(bad+1:end) = [];
      U10r(bad+1:end) = [];
      Tbulk(bad+1:end) = [];
      Rain(bad+1:end) = [];
      Tme(bad+1:end) = [];


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
   nmnTair(i) = nmean(Bins(i).Tair);
   nsdTair(i) = nstd(Bins(i).Tair);
   nmnTbulk(i) = nmean(Bins(i).Tbulk);
   nsdTbulk(i) = nstd(Bins(i).Tbulk);
   nmnU10r(i) = nmean(Bins(i).U10r);
   nsdU10r(i) = nstd(Bins(i).U10r);
   nmnSST(i) = nmean(Bins(i).SST);
   nsdSST(i) = nstd(Bins(i).SST);
end







figure
orient tall
ax(1) = subplot(4,1,1)
plot(bins,nmnTair,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTair(i)-nsdTair(i) nmnTair(i)+nsdTair(i)],'k');
end
title('Dynamo Rain Events - Events > 5 Hours Only','FontSize',16)
ylabel('\DeltaT_{air} [\circC]','FontSize',14)
set(gca,'FontSize',12)
ylim([-3 1])
hline(0)


ax(2) = subplot(4,1,2)
plot(bins,nmnU10r,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnU10r(i)-nsdU10r(i) nmnU10r(i)+nsdU10r(i)],'k');
end
ylabel('\DeltaU_{r_{10}} [m/s]','FontSize',14)
set(gca,'FontSize',12)
ylim([-5 8])
hline(0)

ax(3) = subplot(4,1,3)
plot(bins,nmnTbulk,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTbulk(i)-nsdTbulk(i) nmnTbulk(i)+nsdTbulk(i)],'k');
end
ylabel('\DeltaT_{bulk} [\circC]','FontSize',14)
set(gca,'FontSize',12)
ylim([-0.5 0.2])
hline(0)


ax(4) = subplot(4,1,4)
plot(bins,nmnSST,'k-o','MarkerSize',5','MarkerFaceColor','r','MarkerEdgeColor','r','linewidth',1)
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnSST(i)-nsdSST(i) nmnSST(i)+nsdSST(i)],'k');
end
ylabel('\DeltaSST_{skin} [\circC]','FontSize',14)
xlabel('(T-T_0)/t_{*SST}','FontSize',14)
ylim([-0.8 0.4])
hline(0)
set(gca,'FontSize',12,'ytick',[-0.8:0.4:0.4])
xlim([-2 4])
linkaxes(ax,'x')

print -dpng Dynamo_RainEvent_Composite_long.png
print -dpsc Dynamo_RainEvent_Composite_long.ps


CL(1,:) = [0 0 1];
CL(2,:) = [1 0 0];
CL(3,:) = [0 0 0];
CL(4,:) = [0.5 0 1];

figure
orient tall

ax(1) = subplot(4,1,1)
for j = 1:length(long_events)
   i = long_events(j);
%   plot(Vars(i).TRel11*24*60./Vars(i).Time_to_Max, Vars(i).Tair,'color',CL(j,:))
   plot(Vars(i).TRel11*24*60, Vars(i).Tair-Vars(i).Tair(1+beg_pad1),'color',CL(j,:))
   hold on
end
xlabel('(t-t_0) [minutes]','FontSize',14)
ylabel('\DeltaT_{air} [\circC]','FontSize',14)
hline(0)
title('Dynamo Rain Events > 5 Hours','FontSize',16)


ax(2) = subplot(4,1,2)

for j = 1:length(long_events)
   i = long_events(j);
%   plot(Vars(i).TRel11*24*60./Vars(i).Time_to_Max, Vars(i).Wind,'color',CL(j,:))
   plot(Vars(i).TRel11*24*60, Vars(i).Wind-Vars(i).Wind(1+beg_pad1),'color',CL(j,:))
   hold on
end
xlabel('(t-t_0) [minutes]','FontSize',14)
ylabel('\DeltaU_{10r} [m/s]','FontSize',14)
hline(0)


ax(3) = subplot(4,1,3)
for j = 1:length(long_events)
   i = long_events(j);
%   plot(Vars(i).TRel11*24*60./Vars(i).Time_to_Max, Vars(i).Tbulk,'color',CL(j,:))
   plot(Vars(i).TRel11*24*60, Vars(i).Tbulk-Vars(i).Tbulk(1+beg_pad1),'color',CL(j,:))
   hold on
end
xlabel('(t-t_0) [minutes]','FontSize',14)
ylabel('\DeltaT_{bulk} [\circC]','FontSize',14)
hline(0)


ax(4) = subplot(4,1,4)
indx1 = 1+beg_pad2:-1:1+beg_pad2-18;
for j = 1:length(long_events)
   i = long_events(j);
   SST = Vars(i).SST_smoothed;
   SST = SST - nmean(SST(indx1));
   plot(Vars(i).TRel21*24*60,SST,'color',CL(j,:))
   hold on
end
hline(0)
ylabel('\DeltaSST_{skin} [\circC]','FontSize',14)
xlabel('(t-t_0) [minutes]','FontSize',14)
linkaxes(ax,'x')
xlim([-50 950])


print -dpng Dynamo_RainEvents_LongEvents.png
print -dpsc Dynamo_RainEvents_LongEvents.ps

