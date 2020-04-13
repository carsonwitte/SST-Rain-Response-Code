

% -----------------------------------------------
% Modified to save heat fluxes per Zappa
% LeBel 06.04.2013
%
% -----------------------------------------------

load SST_rain_comparison_peak_data.mat

loc = '/local/data/deadshot1/Analysis/DYNAMO/KT15/IRSST/';
do_compute = 1;



if do_compute

   npeaks = size(peakpos,1);
   %beg_pad = 1*6*10;       % 10 minutes
   beg_pad = 1*6*20;       % 20 minutes
   %end_pad = 1*6*15;	 % 15 minutes
   %end_pad = 1*6*45;	 % 45 minutes
   end_pad = 1*6*60;	 % 60 minutes
   beg_pad1 = 2;
   end_pad1 = 4;

   for i = 1:npeaks
      strtT = peakpos(i,2);	% time
      pkT = peakpos(i,4);		% position
      ndT = peakpos(i,6);		% time
      mm = find(yday >= yday(strtT) & yday <= yday(ndT));		% Edson's data - 10-min avg
      mn = find(yearday >= yday(strtT) & yearday <= yday(ndT));		% Skin SST - 10-sec avg.
      if ~isempty(mn)
       nn = [mn(1)-beg_pad-1:mn(1)-1];
       nn = [nn mn];
       nn = [nn mn(end)+1:mn(end)+1+end_pad];
      else
         nn = mn;
      end
      MM = mm(1)-beg_pad1:mm(end)+end_pad1;
  
      Vars(i).Time_of_Peak = yday(pkT);
      Vars(i).Peak = P(pkT);
      Vars(i).time_indices = MM;
      Vars(i).time_indices_desc = 'position in vector';
      Vars(i).Time = yday(MM); 
      Vars(i).Tair = T10(MM);
      Vars(i).Tbulk = TSea(MM);
      Vars(i).Wind = Ur10(MM);
      Vars(i).Rain = P(MM);
      Vars(i).Rain_desc = 'mm/hr';
      Vars(i).Total_Rain = nsum(P(mm));
      Vars(i).LHF = LHF(MM);
      Vars(i).LHF_desc = 'latent heat flux';
      Vars(i).SHF = SHF(MM);
      Vars(i).SHF_desc = 'sensible heat flux';
      Vars(i).RHF = RHF(MM);
      Vars(i).RHF_desc = 'latent heat flux due to rain';
      Vars(i).OLR = OLR(MM);
      Vars(i).OLR_desc = 'IRup + IRdn';
      Vars(i).Solar = nmean(Solar(MM)); 
      if ~isempty(nn)
         Vars(i).SST = SkinSST(nn);
         Vars(i).SST_desc = 'From KT15';
         Vars(i).SST_smoothed = moving_average(moving_average(moving_average(SkinSST(nn),3),3),3);
         Vars(i).SST_smoothed_desc = 'smoothed version for plotting and analysis';
         Vars(i).time_indices2 = nn;
         Vars(i).time_indices2_desc = 'position in vector';
         Vars(i).T2 = yearday(nn);
         Vars(i).T2_desc = 'Times for SST';
         Vars(i).TRel1 = yday(MM)-yday(pkT);
         Vars(i).TRel1_desc = 'Time relative to peak for Edson''s 10-min data';
         Vars(i).TRel11 = yday(MM)-yday(mm(1));
         Vars(i).TRel11_desc = 'Time relative to beginning of event for Edson''s 10-min data';

         Vars(i).TRel2 = yearday(nn)-( yday(pkT) );
         Vars(i).TRel2_desc = 'Time relative to peak for Zappa''s 10-sec data';
         Vars(i).TRel21 = yearday(nn)-( yday(mm(1)) );
         Vars(i).TRel21_desc = 'Time relative to beginning for Zappa''s 10-sec data';
      else
         Vars(i).SST = NaN;
         Vars(i).SST_desc = 'From KT15';
         Vars(i).SST_smoothed = NaN;
         Vars(i).SST_smoothed_desc = 'smoothed version for plotting and analysis';
         Vars(i).time_indices2 = NaN;
         Vars(i).time_indices2_desc = 'position in vector';
         Vars(i).T2 = NaN;
         Vars(i).T2_desc = 'Times for SST';
         Vars(i).TRel1 = yday(MM)-yday(pkT);
         Vars(i).TRel1_desc = 'Time relative to peak for Edson''s 10-min data';
         Vars(i).TRel11 = yday(MM)-yday(MM(1));
         Vars(i).TRel11_desc = 'Time relative to beginning for Edson''s 10-min data';
         Vars(i).TRel2 = NaN;
         Vars(i).TRel2_desc = 'Time relative to peak for Zappa''s 10-sec data';
         Vars(i).TRel21 = NaN;
         Vars(i).TRel21_desc = 'Time relative to beginning for Zappa''s 10-sec data';
      end
   end


   save SST_rain_comparison_variables_60minpad.mat

else

   load SST_rain_comparison_variables.mat
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


   ax(1) = axes('position',[0.13 0.739 0.775 0.09])
   box on
   mx = -32728;
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      hold on
      Y = Vars(Vindx).Wind;
      Y = Y-Y(1+beg_pad1);
      X = Vars(Vindx).TRel1*24*60;
      %plot(Vars(Vindx).TRel1*24*60, Y,'.-','color',CC(j,:) );
      plot(X, Y,'.-','color',CC(j,:) );
      hold on
      %plot( Vars(Vindx).TRel1(end-end_pad1)*24*60, Y(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      %plot( Vars(Vindx).TRel1(1+beg_pad1)*24*60, Y(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
 
      plot( X(end-end_pad1),Y(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot( X(1+beg_pad1), Y(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )

      %mx1 = ceil(max(Vars(Vindx).Wind)/5);
      mx1 = ceil(max(Y));
      if mx1 > mx
         mx = mx1;
      end
   end
   vline(0)
   hline(0)
   ylabel('Wind Speed [m/s]','FontSize',14)
   ylim([-mx mx])
   title(['Dynamo Rain Event ' num2str( (pg-1)*5 + 1) ' - ' num2str( (pg-1)*5 + 5)],'FontSize',14)


   %ax(2) = subplot(5,1,2)
   ax(2) = axes('position',[0.13 0.611 0.775 0.09])
   box on
   mx = 32728;
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Y = Vars(Vindx).Tair-Vars(Vindx).Tair(beg_pad1+1);
      X = Vars(Vindx).TRel1*24*60;
      %plot(Vars(Vindx).TRel1*24*60, Y,'.-','color',CC(j,:) )
      plot(X, Y,'.-','color',CC(j,:) )
      hold on
      %plot(Vars(Vindx).TRel1(end-end_pad1)*24*60, Y(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      %plot(Vars(Vindx).TRel1(1+beg_pad1)*24*60, Y(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot(X(end-end_pad1), Y(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot(X(1+beg_pad1), Y(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      if min(Y) < mx
         mx = min(Y);
      end
   end
   mx = floor(mx);
   if mod(mx,2) == 1
      mx = mx-1;
   end
   ylim([mx -mx])
   vline(0)
   %vline( (Vars(Vindx).TRel1(1+beg_pad1))*24*60)
   %vline( (Vars(Vindx).TRel1(end-end_pad1))*24*60)
   hline(0)
   ylabel('T_{air} [\circC]','FontSize',14)



   %ax(3) = subplot(5,1,3)
   ax(3) = axes('position',[0.13 0.430 0.775 0.13])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      X = Vars(Vindx).TRel1*24*60;
      hold on
      %plot(Vars(Vindx).TRel1*24*60, Vars(Vindx).Rain,'.-','color',CC(j,:) )
      %plot(Vars(Vindx).TRel1(end-end_pad1)*24*60, Vars(Vindx).Rain(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      %plot(Vars(Vindx).TRel1(1+beg_pad1)*24*60, Vars(Vindx).Rain(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot(X, Vars(Vindx).Rain,'.-','color',CC(j,:) )
      plot(X(end-end_pad1), Vars(Vindx).Rain(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot(X(1+beg_pad1), Vars(Vindx).Rain(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      mx1 = ceil(max(Vars(Vindx).Rain)/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([-2.5 10*mx])
   vline(0)
   %vline( (Vars(Vindx).TRel1(1+beg_pad1))*24*60)
   %vline( (Vars(Vindx).TRel1(end-end_pad1))*24*60)
   ylabel('Rain Rate [mm/hr]','FontSize',14)


%ax(4) = subplot(5,1,4)
   ax(4) = axes('position',[0.13 0.281 0.775 0.10])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      %hold on
      %plot(Vars(Vindx).TRel1(2:end-1), cumsum(Vars(Vindx).Rain(2:end-1))*1/6*(size(Vars(Vindx).Rain,2)-2),'.-','color',CC(j,:) )
      semilogy(Vars(Vindx).TRel1(beg_pad1+1:end-end_pad1)*24*60, cumsum(Vars(Vindx).Rain(beg_pad1+1:end-end_pad1))*1/6*(size(Vars(Vindx).Rain,2)-(beg_pad1+end_pad1)),'.-','color',CC(j,:) )
      hold on
      mx1 = ceil((nsum(Vars(Vindx).Rain)*1/6*(size(Vars(Vindx).Rain,2)-(beg_pad1+end_pad1)))/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([1e-2 mx]*10)
   vline(0)
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
         X = Vars(Vindx).TRel2*24*60;
         %plot(Vars(Vindx).TRel2*24*60, SSTtmp,'color',CC(j,:) )
         %plot(Vars(Vindx).TRel2(end-end_pad)*24*60, SSTtmp(end-end_pad),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
         %plot(Vars(Vindx).TRel2(1+beg_pad)*24*60, SSTtmp(1+beg_pad),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
         plot(X, SSTtmp,'color',CC(j,:) )
         plot(X(end-end_pad), SSTtmp(end-end_pad),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
         plot(X(1+beg_pad), SSTtmp(1+beg_pad),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
      end
      mny = floor(min(SSTtmp)/0.5);
      mxy = ceil(max(SSTtmp)/0.5);
      mxx = max(Vars(Vindx).TRel2(end)*24*60);
      mnx = min(Vars(Vindx).TRel2(1)*24*60);
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
   %XL(2) = ceil(( MxX -1)*10000)/10000+1;
   %XL(1) = floor(( MnX -1)*10000)/10000+1;
   XL(2) = ceil( MxX/10)*10;
   XL(1) = floor( MnX/10)*10;

   ylabel('SST [\circC]','FontSize',14)
   xlabel('Time Relative to Peak Rain Rate [minutes]','FontSize',14)
   ylim([MnY MxY]*0.5)
   vline(0)
   %vline( (Vars(Vindx).TRel2(end-end_pad))*24*60)
   %vline( (Vars(Vindx).TRel2(1+beg_pad))*24*60)
   hline(0)
   set(ax,'xlim',XL); 
   %linkaxes(ax,'x')
   
   outfile = ['Dynamo_Rain_Event_' num2str( (pg-1)*5 + 1) '-' num2str( (pg-1)*5 + 5) '_relpeak_RTUSST'];
   %eval(['print -dpng ' outfile '.png'])
   %eval(['print -dpsc ' outfile '.ps'])
end




return



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
      Tstar = Vars(Vindx).TRel1+Vars(Vindx).Time_of_Peak;
      Tstar = Tstar - Tstar(1+beg_pad1);
      hold on
      Y = Vars(Vindx).Wind-Vars(Vindx).Wind(beg_pad1+1);
      plot(Tstar*24*60, Y,'.-','color',CC(j,:) )
      plot(Tstar(end-end_pad1)*24*60, Y(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot(Tstar(1+beg_pad1)*24*60, Y(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      %mx1 = ceil(max(Vars(Vindx).Wind)/5);
      mx1 = ceil(max(Y));
      if mx1 > mx
         mx = mx1;
      end
   end
   vline(0)
   hline(0)
   ylabel('Wind Speed [m/s]','FontSize',14)
   ylim([-mx mx])
   title(['Dynamo Rain Event ' num2str( (pg-1)*5 + 1) ' - ' num2str( (pg-1)*5 + 5)],'FontSize',14)


   %ax(2) = subplot(5,1,2)
   ax(2) = axes('position',[0.13 0.611 0.775 0.09])
   box on
   mx = 32278;
   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1+Vars(Vindx).Time_of_Peak;
      Tstar = Tstar - Tstar(1+beg_pad1);
      hold on
      Y = Vars(Vindx).Tair-Vars(Vindx).Tair(1+beg_pad1);
      plot(Tstar*24*60, Y, '.-','color',CC(j,:) )
      plot(Tstar(end-end_pad1)*24*60, Y(end-end_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot(Tstar(1+beg_pad1)*24*60, Y(1+beg_pad1), 'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      if min(Y) < mx
         mx = min(Y);
      end
   end
%   ylim([22 30])
   mx = floor(mx);
   if mod(mx,2) == 1
      mx = mx-1;
   end
   ylim([mx -mx])
   vline(0)
   hline(0)
   ylabel('T_{air} [\circC]','FontSize',14)




   %ax(3) = subplot(5,1,3)
   ax(3) = axes('position',[0.13 0.430 0.775 0.13])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1;
      Tstar = Tstar + Vars(Vindx).Time_of_Peak;
      Tstar = Tstar - Tstar(1+beg_pad1);
      hold on
      plot(Tstar*24*60, Vars(Vindx).Rain,'.-','color',CC(j,:) )
      plot(Tstar(end-end_pad1)*24*60, Vars(Vindx).Rain(end-end_pad1) ,'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      plot(Tstar(1+beg_pad1)*24*60, Vars(Vindx).Rain(1+beg_pad1) ,'o','MarkerSize',5, 'MarkerEdgeColor',CC(j,:) )
      mx1 = ceil(max(Vars(Vindx).Rain)/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([-2.5 10*mx])
   vline(0)
   ylabel('Rain Rate [mm/hr]','FontSize',14)


%ax(4) = subplot(5,1,4)
   ax(4) = axes('position',[0.13 0.281 0.775 0.10])
   box on
   mx = -32728;

   for j = 1:nev
      Vindx = (pg-1)*5 + j;
      Tstar = Vars(Vindx).TRel1+Vars(Vindx).Time_of_Peak;
      Tstar = Tstar - Tstar(beg_pad1+1);
      %%plot(Tstar(beg_pad1+1:end-end_pad1), cumsum(Vars(Vindx).Rain(beg_pad1+1:end-end_pad1))*1/6*(size(Vars(Vindx).Rain,2)-2),'.-','color',CC(j,:) )
      %semilogy(Tstar(beg_pad1+1:end-1)*24*60, cumsum(Vars(Vindx).Rain(beg_pad1+1:end-end_pad1))*1/6*(size(Vars(Vindx).Rain,2)-2),'.-','color',CC(j,:) )
      semilogy(Tstar(beg_pad1+1:end-end_pad1)*24*60, cumsum(Vars(Vindx).Rain(beg_pad1+1:end-end_pad1))*1/6*(size(Vars(Vindx).Rain,2)-(beg_pad1+end_pad1)),'.-','color',CC(j,:) )
      hold on
 
      mx1 = ceil((nsum(Vars(Vindx).Rain)*1/6*(size(Vars(Vindx).Rain,2)-2))/10);
      if mx1 > mx
         mx = mx1;
      end
   end
   ylim([1e-2 mx]*10)
   vline(0)
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
         Tstar = Vars(Vindx).TRel2+Vars(Vindx).Time_of_Peak;
         Tstar = Tstar - Tstar(1+beg_pad);
         SSTtmp = SSTtmp-nmean(SSTtmp(1*6*9:1*6*10));
         hold on
         plot(Tstar*24*60, SSTtmp,'color',CC(j,:) )
         plot(Tstar(end-end_pad)*24*60, SSTtmp(end-end_pad),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
         plot(Tstar(1+beg_pad)*24*60, SSTtmp(1+beg_pad),'o','MarkerSize',4,'MarkerFaceColor',CC(j,:),'MarkerEdgeColor',CC(j,:) )
      end

      mny = floor(min(SSTtmp)/0.5);
      mxy = ceil(max(SSTtmp)/0.5);
      %mxx = max(Vars(Vindx).TRel2(end));
      %mnx = min(Vars(Vindx).TRel2(1));
      mxx = max(Tstar);
      mnx = min(Tstar);
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
   %XL(2) = ceil(( MxX -1)*10000)/10000+1;
   %%XL(1) = floor(( MnX -1)*10000)/10000+1;
   %XL(1) = 0.9999

   ylabel('SST [\circC]','FontSize',14)
   xlabel('Time Relative to Beginning of Event [minutes]','FontSize',14)
   ylim([MnY MxY]*0.5)
   vline(0)
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




