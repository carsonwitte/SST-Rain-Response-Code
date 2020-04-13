%
% Basic rainfall peak detection and collation of Edson's met data
% Reversed rainfall series to detect the end of each event as well.
% 
% First step in the process!


loc = 'C:\Users\ZappaLab\Google Drive\Experiments\2020 SST Rain Response Paper\Data\Working';

load([loc '\KT15_corrected_SST_10sec_for_distribution.mat'])

SkinSST = SST;
yearday = yday;

load ../Revelle10minutesLeg2_r1.mat
nn = find(P < 1e-6);
P(nn) = 0;

L1 = 4249;
L2 = 4321;
L3 = 2161;


ydaytot = [];
Ur10tot = [];
Ptot = [];
T10tot = [];
RHFtot = [];
SHFtot = [];
LHFtot = [];
IRtot = [];
TSeatot = [];
Solartot = [];

ydaytot = [ydaytot yday];
Ur10tot = [Ur10tot Ur10];
T10tot = [T10tot T10];
Ptot = [Ptot P];
RHFtot = [RHFtot rhf];
LHFtot = [LHFtot lhf];
SHFtot = [SHFtot shf];
IRtot = [IRtot IRup+IRdn];
TSeatot = [TSeatot Tsea];
Solartot = [Solartot Solardn];


% find rainfall peaks
thresh = 0.45;			% mm/hr, change value to define peak

[maxval,minval,minval2] = peakdet2(P,0.45);	% minval gives position & value of first zero following peak;
peakpos =  ones(size(maxval,1),7);
size(peakpos)
MM = size(maxval,1);

% value of peak
peakpos(1:MM,1) = maxval(:,2);
size(peakpos)

% position of peak
peakpos(1:MM,4) = maxval(:,1);

% time of peak
peakpos(1:MM,5) = yday(maxval(:,1));

%following zero
peakpos(1:MM-1,6) = minval(:,1);
peakpos(MM,6) = NaN;

peakpos(1:MM-1,7) = yday(minval(:,1));
peakpos(MM,7) = NaN;

%leading zero
peakpos(2:MM,2) = sort(minval2(:,1));
peakpos(1,2) = NaN;

peakpos(2:MM,3) = yday(sort(minval2(:,1)));
peakpos(1,3) = NaN;

% hand-edit values
%
% first zero for Leg II

peakpos(1,2) = 46;
peakpos(1,3) = yday(46);

peakpos(59,6) = 4190;
peakpos(59,7) = yday(4190);




for j = 2:3
   if j == 2
      Padd = L1;
   else
      Padd = L1 + L2;
   end

   M0 = size(peakpos,1);
   file = ['../Revelle10minutesLeg' num2str(j+1) '_r1.mat']
   load(file)
   ydaytot = [ydaytot yday];
   Ur10tot = [Ur10tot Ur10];
   T10tot = [T10tot T10];
   Ptot = [Ptot P];
   TSeatot = [TSeatot Tsea];
   RHFtot = [RHFtot rhf];
   LHFtot = [LHFtot lhf];
   IRtot = [IRtot IRup+IRdn];
   SHFtot = [SHFtot shf];
   Solartot = [Solartot Solardn];
   clear maxval minval Z ZZ
   nn = find(P < 1e-6);
   P(nn) = 0;

   [maxval,minval,minval2] = peakdet2(P,0.45);	% minval gives position & value of first zero following peak;
   MM = size(maxval,1);
   if j == 2
      minval(MM,1) = 4320;
      minval(MM,2) = P(4320);
      minval2(MM,1) = 6;
      minval2(MM,2) = P(6);
   else
      minval(MM,1) = 1185;
      minval(MM,2) = P(1185);
      minval2(MM,1) = 124;
      minval2(MM,2) = P(124);
   end

 
   peakpos(M0+1:M0+MM,1) = maxval(:,2);		% peak value
   peakpos(M0+1:M0+MM,4) = maxval(:,1)+Padd;	% peak position
   peakpos(M0+1:M0+MM,5) = yday(maxval(:,1));	% peak position - yearday


   %following zero
   peakpos(M0+1:M0+MM,6) = minval(:,1)+Padd;
   peakpos(M0+1:M0+MM,7) = yday(minval(:,1));

   %following zero
   peakpos(M0+1:M0+MM,2) = sort(minval2(:,1))+Padd;
   peakpos(M0+1:M0+MM,3) = yday(sort(minval2(:,1)));


end


clear wdir wdirR rhf shf lhf nn yday j i interped Precip P10 SOG COG 
clear SSQ SST Solarup Solardn IRup IRdn Heading Evap E RH10 Qsea Q10 Ur10 cdir cspd 
clear Lat Lon MM M0 Padd Pair10 T10 U10 stress Tsea IRup IRdn 

T10 = T10tot;
Ur10 = Ur10tot;
P = Ptot;
yday = ydaytot;
Solar = Solartot;
RHF = RHFtot;
SHF = SHFtot;
LHF = LHFtot;
OLR = IRtot;
TSea = TSeatot;
clear T10tot Ur10tot Ptot ydaytot TSeatot RHFtot Solartot TSeatot LHFtot SHFtot IRtot



save SST_rain_comparison_peak_data.mat



peakpos_all = peakpos;

tmp = (peakpos(1:end-1,6)-peakpos(2:end,2));
dupes = find(tmp == 0);				% end of last peak is beginning of next peak
ht1 = P(peakpos(dupes,4) );
ht2 = P(peakpos(dupes+1,4) );


while(length(dupes) > 0)
   length(dupes)
   [dum,minpk] = min([ht1(1) ht2(1)]);
   [dum,maxpk] = max([ht1(1) ht2(1)]);
   rowrm = dupes(1) -1 + minpk; 
   rowkeep = dupes(1) -1 + maxpk; 
   %peakpos(row,:) = [];
   if rowkeep > rowrm			% keep second peak 
      tmprow = peakpos(rowkeep,:);
      tmprow(:,2:3) = peakpos(rowrm,2:3);		% save beginning
      peakpos(rowkeep,:) = tmprow;
   else					% keep first peak
      tmprow = peakpos(rowkeep,:);
      tmprow(:,6:7) = peakpos(rowrm,6:7);   % save ending
   end
   peakpos(rowkeep,:) = tmprow;
   peakpos(rowrm,:) = [];
   
   tmp = (peakpos(1:end-1,6)-peakpos(2:end,2));
   dupes = find(tmp == 0);				% end of last peak is beginning of next peak
   ht1 = P(peakpos(dupes,4) );
   ht2 = P(peakpos(dupes+1,4) );
end



clear tmp ht1 ht2 rowkeep rowrm dum minpk maxpk tmprow row 

save SST_rain_comparison_peak_data.mat




