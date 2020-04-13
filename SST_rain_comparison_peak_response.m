
%  Calculate time of peak SST response 
%  A bit of a hack still
%  Consider doing piece-wise fit of data to a line
%
%

load SST_rain_comparison_peak_data.mat

loc = '/local/data/deadshot1/Analysis/DYNAMO/KT15/IRSST/';
%load SST_rain_comparison_variables.mat
load SST_rain_comparison_variables_60minpad.mat

npeaks = size(Vars,2);
   beg_pad2 = 1*6*20;    % 20 minutes	SST
%   end_pad2 = 1*6*45;	 % 45 minutes	SST
   end_pad2 = 1*6*60;	 % 60 minutes	SST
   beg_pad1 = 2;	 % 20 minutes  Edson
   end_pad1 = 4;	 % 40 minutes  Edson


do_compute = 1;


if do_compute
   for i = 1:npeaks

      %T = Vars(i).TRel2;
      T = Vars(i).T2;
      Tmp = Vars(i).SST;
      Tmp_smoothed = Vars(i).SST_smoothed;
      Tpk = Vars(i).Time_of_Peak;
      Tair = Vars(i).Tair;


      % time of actual rain event				% just look for response within rain event
      mn = find(Vars(i).T2 >= Vars(i).Time(beg_pad1+1) & Vars(i).T2 <= Vars(i).Time(end-end_pad1) );
      mn2 = find(Vars(i).T2 >= Vars(i).Time(beg_pad1+2) & Vars(i).T2 <= Vars(i).Time(end-end_pad1) );

      if ~isempty(mn)
         SST0 = nmean(Tmp_smoothed(mn(1)-12:mn(1)+5) );	   % 3 minutes before rain event
         %Tnot = T(mn(1));					% beginning of rain event
         Tnot = T(mn2(1));					% beginning of rain event
         %mx = find(SST0 - Tmp_smoothed(mn) == max(SST0 - Tmp_smoothed(mn)) );
         mx = find(SST0 - Tmp_smoothed(mn2) == max(SST0 - Tmp_smoothed(mn2)) );
         %Tmax = T(mn(mx)) ; 				   % time of max response
         Tmax = T(mn2(mx)) ; 				   % time of max response
         %Ttot = (T(mn(end)) - T(mn(1)))*24*60;		   % length of rain event
         Ttot = (T(mn2(end)) - T(mn2(1)))*24*60; 	   % length of rain event
         Tot = (Vars(i).Time(end-end_pad1) - Vars(i).Time(beg_pad1+2))*24*60;
         % find time of recovery 
         nn = find(T > Tmax);                              % times after SST response
         bb = find(Vars(i).SST_smoothed(nn) >= SST0);
         if ~isempty(bb)
            indx = nn(bb(1));
         end
         

         %fprintf('%5d %6.2f %6.2f %6.2f %6.2f\n', i, Ttot,(Tpk-Tnot)*24*60, (Tmax-Tnot)*24*60, (Tmax - Tpk)*24*60 );
         Vars(i).SST0 = SST0;
         Vars(i).SST0_desc = '3 minute average just prior to rain event';
         Vars(i).Max_Response = Tmp_smoothed(mn2(mx));
         Vars(i).Max_Response_desc = 'minimum SST during rain event';
         jkl = find(Tair(1+beg_pad1:end-end_pad1) == min(Tair(1+beg_pad1:end-end_pad1)) );
         Vars(i).Max_Response_Tair = Tair(beg_pad1+jkl);
         Vars(i).Max_Response_Tair_desc = 'minimum air temperature during rain event';
         Vars(i).Tair_response =  Tair(beg_pad1+jkl) -Tair(beg_pad1+1);
         Vars(i).Tair_response_desc = 'Change in air temperature (degrees C)';
         Vars(i).Time_of_Max_Response = Tmax;
         Vars(i).Time_to_Max = (Tmax-Tnot)*24*60;
         Vars(i).Time_to_Max_desc = 'Time to reach max SST response';
         Vars(i).Time_between_peak_max = (Tpk-Tmax)*24*60;
         Vars(i).Time_between_peak_max_desc = 'minutes between max response and peak rain rate';
         Vars(i).Length_of_Rain_Event = Ttot;
         Vars(i).Length_of_Rain_Event_desc = 'minutes';
         if ~isempty(bb)
            Vars(i).Time_of_Recovery = Vars(i).T2(indx);
            Vars(i).Time_of_Recovery_desc = 'Time of return of SST to SST0';
            Vars(i).Time_to_Recovery = (Vars(i).T2(indx) - Tmax)*24*60;
            Vars(i).Time_to_Recovery_desc = 'Minutes between maximum SST response and return to SST0';
         else
            Vars(i).Time_of_Recovery = NaN;
            Vars(i).Time_of_Recovery_desc = 'Time of return of SST to SST0';
            Vars(i).Time_to_Recovery = NaN;
            Vars(i).Time_to_Recovery_desc = 'Minutes between maximum SST response and return to SST0';
         end
      else
         Vars(i).SST0 = NaN;
         Vars(i).SST0_desc = NaN;
         Vars(i).Max_Response = NaN;
         Vars(i).Max_Response_desc = NaN;
         Vars(i).Time_of_Max_Response = NaN;
         Vars(i).Time_to_Max = NaN;
         Vars(i).Time_to_Max_desc = NaN;
         Vars(i).Time_between_peak_max = NaN;
         Vars(i).Time_between_peak_max_desc = NaN;
         Vars(i).Length_of_Rain_Event = NaN;
         Vars(i).Length_of_Rain_Event_desc = NaN;
         Vars(i).Time_of_Recovery = NaN;
         Vars(i).Time_of_Recovery_desc = 'Time of return of SST to SST0';
         Vars(i).Time_to_Recovery = NaN;
         Vars(i).Time_to_Recovery_desc = 'Minutes between maximum SST response and return to SST0';
      end
   end


   save SST_rain_comparison_variables2_60minpad.mat Vars beg_pad2 end_pad2 beg_pad1 end_pad1

else		% do_compute = 0

   load SST_rain_comparison_variables2_60minpad.mat
end



return

nev = 5;
npages = floor(npeaks/nev);
CC(1,:) = [0 0 1];
CC(2,:) = [1 0 0];
CC(3,:) = [0 1 0];
CC(4,:) = [0 0 0];
CC(5,:) = [0.7 0 1];

for pg = 1:npages
   figure
   orient landscape 


   ax(1) = subplot(3,1,1)
   box on
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1;
      Tstar = Tstar * Vars(Vindx).Time_of_Peak;
      Tstar = Tstar / Tstar(2);
      hold on
      plot(Tstar, Vars(Vindx).Tair-Vars(Vindx).Tair(1),'.-','color',CC(j,:) )
   end
%   ylim([22 30])
   vline(1)
   hline(0)
   ylabel('T_{air} [\circC]','FontSize',14)
   title(['Dynamo Rain Event ' num2str( (pg-1)*5 + 1) ' - ' num2str( (pg-1)*5 + 5)],'FontSize',14)


   ax(2) = subplot(3,1,2)
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1;
      Tstar = Tstar * Vars(Vindx).Time_of_Peak;
      Tstar = Tstar / Tstar(2);
      hold on
      plot(Tstar, Vars(Vindx).Rain,'.-','color',CC(j,:) )
      mx1 = ceil(max(Vars(Vindx).Rain)/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([0 mx]*10)
   vline(1)
   ylabel('Rain Rate [mm/hr]','FontSize',14)



   ax(3) = subplot(3,1,3)
   box on
   MnY = 32728;
   MxY = -32728;
   MnX = 32728;
   MxX = -32728;
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      %SSTtmp = Vars(Vindx).SST;
      SSTtmp = Vars(Vindx).SST_smoothed;
      if ~isnan(SSTtmp)
         Tstar = Vars(Vindx).TRel2;
         Tstar = Tstar * Vars(Vindx).Time_of_Peak;
         %Tstar = Tstar / Tstar(1*6*10+1);
         SSTtmp = SSTtmp-nmean(SSTtmp(1*6*9:1*6*10));
         hold on
         plot(Tstar/Tstar(1*6*10+1), SSTtmp,'color',CC(j,:) )
         plot(Tstar(end-1*6*15)/Tstar(1*6*10+1), SSTtmp(end-1*6*15),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
         vline(Vars(Vindx).Time_of_Max_Response/Tstar(1*6*10+1),CC(j,:) );
      end

      mny = floor(min(SSTtmp)/0.5);
      mxy = ceil(max(SSTtmp)/0.5);
      mxx = max(Vars(Vindx).TRel2(end));
      mnx = min(Vars(Vindx).TRel2(1));
      if mny < MnY
         MnY = mny;
      end
      if mxy > MxY
         MxY = mxy;
      end
      if mnx < MnX
         MnX = mnx;
      end
      if mxx > MxX
         MxX = mxx;
      end
   end
   XL(2) = ceil(( MxX -1)*10000)/10000+1;
   %XL(1) = floor(( MnX -1)*10000)/10000+1;
   XL(1) = 0.9999

   ylabel('SST [\circC]','FontSize',14)
   xlabel('Time Relative to Beginning of Event','FontSize',14)
   ylim([MnY MxY]*0.5)
   vline(1)
   hline(0)
   XL = xlim;
   set(ax,'xlim',XL);
   %linkaxes(ax,'x')

   outfile = ['Dynamo_Rain_Event_' num2str( (pg-1)*5 + 1) '-' num2str( (pg-1)*5 + 5) '_with_max_response_RTUSST'];
   eval(['print -dpng ' outfile '.png'])
   eval(['print -dpsc ' outfile '.ps'])
   orient portrait
   set(gcf,'papersize',[11 8.5])
   set(gcf,'paperposition',[.25 .25 10.5 8])
   %eval(['print -dpng ' outfile '.png'])

end



return



%----------------------------------------------------------------
%   Plot individual rain events   - Relative peak rainfall rate
%----------------------------------------------------------------

nev = 5;
npages = floor(npeaks/nev);
CC(1,:) = [0 0 1];
CC(2,:) = [1 0 0];
CC(3,:) = [0 1 0];
CC(4,:) = [0 0 0];
CC(5,:) = [0.7 0 1];

for pg = 1:npages
   figure
   orient tall


   %ax(1) = subplot(5,1,1)
   ax(1) = axes('position',[0.13 0.739 0.775 0.09])
   box on
   mx = -32728;
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      hold on
      plot(Vars(Vindx).TRel1, Vars(Vindx).Wind,'.-','color',CC(j,:) )
      mx1 = ceil(max(Vars(Vindx).Wind)/5);
      if mx1 > mx
         mx = mx1;
      end
   end
   vline(1)
   hline(0)
   ylabel('Wind Speed [m/s]','FontSize',14)
   ylim([0 mx]*5)
   title(['Dynamo Rain Event ' num2str( (pg-1)*5 + 1) ' - ' num2str( (pg-1)*5 + 5)],'FontSize',14)


   %ax(2) = subplot(5,1,2)
   ax(2) = axes('position',[0.13 0.611 0.775 0.09])
   box on
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      hold on
      plot(Vars(Vindx).TRel1, Vars(Vindx).Tair-Vars(Vindx).Tair(1),'.-','color',CC(j,:) )
   end
%   ylim([22 30])
   vline(1)
   hline(0)
   ylabel('T_{air} [\circC]','FontSize',14)



   %ax(3) = subplot(5,1,3)
   ax(3) = axes('position',[0.13 0.430 0.775 0.13])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      hold on
      plot(Vars(Vindx).TRel1, Vars(Vindx).Rain,'.-','color',CC(j,:) )
      mx1 = ceil(max(Vars(Vindx).Rain)/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([0 mx]*10)
   vline(1)
   ylabel('Rain Rate [mm/hr]','FontSize',14)


%ax(4) = subplot(5,1,4)
   ax(4) = axes('position',[0.13 0.281 0.775 0.10])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      %hold on
      %plot(Vars(Vindx).TRel1(2:end-1), cumsum(Vars(Vindx).Rain(2:end-1))*1/6*(size(Vars(Vindx).Rain,2)-2),...
      %        '.-','color',CC(j,:) )
      semilogy(Vars(Vindx).TRel1(2:end-1), cumsum(Vars(Vindx).Rain(2:end-1))*1/6*(size(Vars(Vindx).Rain,2)-2), ...
           '.-','color',CC(j,:) )
      hold on
      mx1 = ceil((nsum(Vars(Vindx).Rain)*1/6*(size(Vars(Vindx).Rain,2)-2))/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([1e-2 mx]*10)
   vline(1)
   ylabel('Volume Rain [mm]','FontSize',14)



%   ax(5) = subplot(5,1,5)
    ax(5) = axes('position',[0.13 0.11 0.775 0.13])
   box on
   MnY = 32728;
   MxY = -32728;
   MnX = 32728;
   MxX = -32728;
  
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      %SSTtmp = Vars(Vindx).SST;
      SSTtmp = Vars(Vindx).SST_smoothed;
      if ~isnan(SSTtmp) 
         SSTtmp = SSTtmp-nmean(SSTtmp(1:18));
         hold on
         plot(Vars(Vindx).TRel2, SSTtmp,'color',CC(j,:) )
         plot(Vars(Vindx).TRel2(end-1*6*15), SSTtmp(end-1*6*15),'o','MarkerSize',3,'MarkerFaceColor',...
             CC(j,:),'MarkerEdgeColor',CC(j,:) )
      end
      mny = floor(min(SSTtmp)/0.5);
      mxy = ceil(max(SSTtmp)/0.5);
      mxx = max(Vars(Vindx).TRel2(end));
      mnx = min(Vars(Vindx).TRel2(1));
      if mny < MnY
         MnY = mny;
      end
      if mxy > MxY
         MxY = mxy;
      end
      if mnx < MnX
         MnX = mnx;
      end
      if mxx > MxX
         MxX = mxx;
      end
   end
   XL(2) = ceil(( MxX -1)*10000)/10000+1;
   XL(1) = floor(( MnX -1)*10000)/10000+1;

   ylabel('SST [\circC]','FontSize',14)
   xlabel('Time Relative to Peak Rain Rate','FontSize',14)
   ylim([MnY MxY]*0.5)
   vline(1)
   hline(0)
   set(ax,'xlim',XL); 
   %linkaxes(ax,'x')
   
   outfile = ['Dynamo_Rain_Event_' num2str( (pg-1)*5 + 1) '-' num2str( (pg-1)*5 + 5) '_relpeak_RTUSST'];
   eval(['print -dpng ' outfile '.png'])
   eval(['print -dpsc ' outfile '.ps'])
end







%-------------------------------------------------------------------
%   Plot individual rain events   - Relative to beginning of event
%-------------------------------------------------------------------

nev = 5;
npages = floor(npeaks/nev);
CC(1,:) = [0 0 1];
CC(2,:) = [1 0 0];
CC(3,:) = [0 1 0];
CC(4,:) = [0 0 0];
CC(5,:) = [0.7 0 1];

for pg = 1:npages
   figure
   orient tall

   %ax(1) = subplot(5,1,1)
   %ax(1) = axes('position',[0.13 0.735 0.775 0.09])
   ax(1) = axes('position',[0.13 0.739 0.775 0.09])
   box on
   mx = -32728;
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1;
      Tstar = Tstar * Vars(Vindx).Time_of_Peak;
      Tstar = Tstar / Tstar(2);
      hold on
      plot(Tstar, Vars(Vindx).Wind,'.-','color',CC(j,:) )
      mx1 = ceil(max(Vars(Vindx).Wind)/5);
      if mx1 > mx
         mx = mx1;
      end
   end
   vline(1)
   hline(0)
   ylabel('Wind Speed [m/s]','FontSize',14)
   ylim([0 mx]*5)
   title(['Dynamo Rain Event ' num2str( (pg-1)*5 + 1) ' - ' num2str( (pg-1)*5 + 5)],'FontSize',14)


   %ax(2) = subplot(5,1,2)
   ax(2) = axes('position',[0.13 0.611 0.775 0.09])
   box on
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1;
      Tstar = Tstar * Vars(Vindx).Time_of_Peak;
      Tstar = Tstar / Tstar(2);
      hold on
      plot(Tstar, Vars(Vindx).Tair-Vars(Vindx).Tair(1),'.-','color',CC(j,:) )
   end
%   ylim([22 30])
   vline(1)
   hline(0)
   ylabel('T_{air} [\circC]','FontSize',14)




   %ax(3) = subplot(5,1,3)
   ax(3) = axes('position',[0.13 0.430 0.775 0.13])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1;
      Tstar = Tstar * Vars(Vindx).Time_of_Peak;
      Tstar = Tstar / Tstar(2);
      hold on
      plot(Tstar, Vars(Vindx).Rain,'.-','color',CC(j,:) )
      mx1 = ceil(max(Vars(Vindx).Rain)/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([0 mx]*10)
   vline(1)
   ylabel('Rain Rate [mm/hr]','FontSize',14)


%ax(4) = subplot(5,1,4)
   ax(4) = axes('position',[0.13 0.281 0.775 0.10])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1(2:end-1);
      Tstar = Tstar * Vars(Vindx).Time_of_Peak;
      Tstar = Tstar / Tstar(2);
      %plot(Tstar, cumsum(Vars(Vindx).Rain(2:end-1))*1/6*(size(Vars(Vindx).Rain,2)-2),'.-','color',CC(j,:) )
      semilogy(Tstar, cumsum(Vars(Vindx).Rain(2:end-1))*1/6*(size(Vars(Vindx).Rain,2)-2),'.-','color',CC(j,:) )
      hold on
 
      mx1 = ceil((nsum(Vars(Vindx).Rain)*1/6*(size(Vars(Vindx).Rain,2)-2))/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([1e-2 mx]*10)
   vline(1)
   ylabel('Volume Rain [mm]','FontSize',14)



%   ax(5) = subplot(5,1,5)
    ax(5) = axes('position',[0.13 0.11 0.775 0.13])
   box on
   MnY = 32728;
   MxY = -32728;
   MnX = 32728;
   MxX = -32728;
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      %SSTtmp = Vars(Vindx).SST;
      SSTtmp = Vars(Vindx).SST_smoothed;
      if ~isnan(SSTtmp)
         Tstar = Vars(Vindx).TRel2;
         Tstar = Tstar * Vars(Vindx).Time_of_Peak;
         Tstar = Tstar / Tstar(1*6*10+1);
         SSTtmp = SSTtmp-nmean(SSTtmp(1*6*9:1*6*10));
         hold on
         plot(Tstar, SSTtmp,'color',CC(j,:) )
         plot(Tstar(end-1*6*15), SSTtmp(end-1*6*15),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
      end

      mny = floor(min(SSTtmp)/0.5);
      mxy = ceil(max(SSTtmp)/0.5);
      mxx = max(Vars(Vindx).TRel2(end));
      mnx = min(Vars(Vindx).TRel2(1));
      if mny < MnY
         MnY = mny;
      end
      if mxy > MxY
         MxY = mxy;
      end
      if mnx < MnX
         MnX = mnx;
      end
      if mxx > MxX
         MxX = mxx;
      end
   end
   XL(2) = ceil(( MxX -1)*10000)/10000+1;
   %XL(1) = floor(( MnX -1)*10000)/10000+1;
   XL(1) = 0.9999

   ylabel('SST [\circC]','FontSize',14)
   xlabel('Time Relative to Beginning of Event','FontSize',14)
   ylim([MnY MxY]*0.5)
   vline(1)
   hline(0)
   %linkaxes(ax,'x')
   XL = xlim;
   set(ax,'xlim',XL); 

   outfile = ['Dynamo_Rain_Event_' num2str( (pg-1)*5 + 1) '-' num2str( (pg-1)*5 + 5) '_RTUSST'];
   eval(['print -dpng ' outfile '.png'])
   eval(['print -dpsc ' outfile '.ps'])
end








return

%-----------------------------------------------------
%   Bin average rain events
%-----------------------------------------------------

maxL = -32768;
minL = 32768;

for i = 1:npeaks
   if length(Vars(i).TRel2) > 0
      if Vars(i).TRel2(end) > maxL
         maxL = Vars(i).TRel2(end);
      end
      if Vars(i).TRel2(1) < minL
         minL = Vars(i).TRel2(1);
      end
   end
end


C = jet(npeaks);



% Do some bin-averaging

bins = 0.9989:0.000025:1.00120;

nmnSST = zeros(length(bins),1);
nsdSST = zeros(length(bins),1);
nmnTair = zeros(length(bins),1);
nsdTair = zeros(length(bins),1);
nmnP = zeros(length(bins),1);
nsdP = zeros(length(bins),1);
nmnU10r = zeros(length(bins),1);
nsdU10r = zeros(length(bins),1);
npnts = zeros(length(bins),1);
npnts2 = zeros(length(bins),1);


% Try using change from initial condition rather than straight value

for j = 1:npeaks
   TRel1 = Vars(j).TRel1;
   TRel2 = Vars(j).TRel2;
   Tair = Vars(j).Tair;
   U10r = Vars(j).Wind;
   SST = Vars(j).SST;
   nn = Vars(j).time_indices2;

   Tair = Tair - Tair(1);
   U10r = U10r - U10r(1);
   if ~isempty(nn)
      SST = SST - nmean(SkinSST(nn(1):nn(1)+10));
   end

   indx2 = NaN*ones(length(Tair),1);
   for k = 1:length(TRel1)
      indx = find(bins > TRel1(k) );
      indx = indx(1)-1;
      indx2(k) = indx;
      nmnTair(indx) = nmnTair(indx) + Tair(k); 
      nmnU10r(indx) = nmnU10r(indx) + U10r(k); 
      npnts(indx) = npnts(indx) + 1;
   end
   Vars(j).indx2 = indx2;
   if ~isempty(nn)
      for k = 1:length(TRel2)
         indx = find(bins > TRel2(k) );
         indx = indx(1)-1;
         nmnSST(indx) = nmnSST(indx) + SST(k); 
         npnts2(indx) = npnts2(indx) + 1;
      end
   end
end
nmnTair = nmnTair./npnts;
nmnU10r = nmnU10r./npnts;
nmnSST = nmnSST./npnts2;


for j = 1:npeaks
   TRel1 = Vars(j).TRel1;
   TRel2 = Vars(j).TRel2;
   Tair = Vars(j).Tair;
   U10r = Vars(j).Wind;
   SST = Vars(j).SST;

   Tair = Tair - Tair(1);
   U10r = U10r - U10r(1);
   if ~isempty(nn)
      SST = SST - nmean(SkinSST(nn(1)-10:nn(1)+10));
   end

   for k = 1:length(TRel1)
      indx = find(bins > TRel1(k) );
      indx = indx(1)-1;
      nsdTair(indx) = nsdTair(indx) + (Tair(k)-nmnTair(indx)).^2; 
      nsdU10r(indx) = nsdU10r(indx) + (U10r(k)-nmnU10r(indx)).^2; 
   end
   if ~isempty(nn)
      for k = 1:length(TRel2)
         indx = find(bins > TRel2(k) );
         indx = indx(1)-1;
         nsdSST(indx) = nsdSST(indx) + (SST(k)-nmnSST(indx)).^2; 
      end
   end
end

nsdTair = nsdTair./(npnts-1);
nsdU10r = nsdU10r./(npnts-1);
nsdSST = nsdSST./(npnts2-1);


figure
subplot(2,1,1)
hold on
for i = 1:npeaks
   plot(Vars(i).TRel1,Vars(i).Rain,'color',C(i,:) );
end




figure
orient tall
ax(1) = subplot(3,1,1)
plot(bins,nmnTair,'-o','MarkerSize',5')
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnTair(i)-nsdTair(i) nmnTair(i)+nsdTair(i)]);
end
title('Dynamo Rain Events','FontSize',14)
ylabel('Tair [\circC]','FontSize',14)


ax(2) = subplot(3,1,2)
plot(bins,nmnU10r,'-o','MarkerSize',5')
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnU10r(i)-nsdU10r(i) nmnU10r(i)+nsdU10r(i)]);
end
ylabel('U_r_{10} [m/s]','FontSize',14)

ax(3) = subplot(3,1,3)
plot(bins,nmnSST,'-o','MarkerSize',5')
hold on
for i = 1:length(bins)
   plot([bins(i) bins(i)],[nmnSST(i)-nsdSST(i) nmnSST(i)+nsdSST(i)]);
end
ylabel('KT15 SST [\circC]','FontSize',14)
xlabel('Relative Time from Peak Rainfall','FontSize',14)
linkaxes(ax,'x')




