function INUdata = readStreamsIMU(filepath,inu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INUDATA = readStreamsIMU(filepath, INU) returns a structure containing
% orientation data from the specified INU device.  It is set up to read raw 
% files with 16-byte header, containing the number of frames and amount of
% offset between the header + timestamps to the actual raw data 
% (should be 16 bytes + NFRAMES bytes, since timestamps are 4 bytes).
%
% INPUT:
%        filepath--string to specify folder where INU raw files are located
%             INU--string to specify 'XSENS', 'MICRO', or 'WATSON' INU 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%         INUdata with fields for XSENS or MICROSTRAIN device:
%             PATH--full path of raw file
%          YMD_UTC--UTC year-month-day-hour-min-sec
%            FRNUM--frame number of corresponding raw data 
%             ROLL--positive rotation clockwise around x-axis,
%                   (wire-to-front), range [-180...+180]
%            PITCH--positive rotation clockwise around y-axis,
%                   (wire-to-left), range [-90...+90]
%              AZM--positive rotation clockwise around z-axis, (wire-to-up),
%                   range [-180...+180]
%            VROLL--roll velocity        (rad/s)
%           VPITCH--pitch veloctiy       (rad/s)
%             VAZM--azimuth velocity     (rad/s)
%            AROLL--roll acceleration    (m/s^2)
%           APITCH--pitch acceleration   (m/s^2)
%             AAZM--azimuth acceleration (m/s^2)
%
%         INUdata with fields for WATSON device:
%             PATH--full path of processed file
%          YMD_UTC--UTC year-month-day-hour-min-sec
%            FRNUM--frame number of corresponding raw data 
%           OMEGA--
%             PHI--
%           KAPPA--
%             VEL--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Scott Brown
% LDEO, Columbia University
% 12 August 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xsens=false; micro=false; watson=false; 
switch upper(inu)
    case 'XSENS'
        xsens=true;
    case 'MICRO'
        micro=true;
    case 'WATSON'
        watson=true;
    otherwise
        error('Enter INU device as ''XSENS'', ''MICRO'', or ''WATSON''.')
end

fid        = fopen(filepath);
nframes    = fread(fid,1,'uint64');
byteOffset = fread(fid,1,'uint64');

%initialize structure fields
INUdata.path   = filepath;
INUdata.YMDHMS = nan(nframes,6);
INUdata.frnum  = (0:nframes-1)';
if xsens || micro
    INUdata.roll    = nan(nframes,1);
    INUdata.pitch   = nan(nframes,1);
    INUdata.azm     = nan(nframes,1);
    INUdata.vroll   = nan(nframes,1);
    INUdata.vpitch  = nan(nframes,1);
    INUdata.vazm    = nan(nframes,1);
    INUdata.aroll   = nan(nframes,1);
    INUdata.apitch  = nan(nframes,1);
    INUdata.aazm    = nan(nframes,1);
elseif watson
    INUdata.omega = nan(nframes,1);
    INUdata.phi   = nan(nframes,1);
    INUdata.kappa = nan(nframes,1);
    INUdata.vel   = nan(nframes,1);
end

%Frame numbers start at 0.
for frIDX = 0:nframes-1
    %seek past 16-byte header to get to timestamps
    fseek(fid,16+frIDX*8,'bof'); 
    tmstmp = datevec(datenum(1601,1,1,0,0,0) +...
             fread(fid,1,'uint64')/(60*60*24*1e7));
    INUdata.YMDHMS(frIDX+1,:) = tmstmp;
    
    if xsens || micro
        %seek past header and timestamps to get to 21-byte data packets
        fseek(fid,byteOffset+21*8*frIDX,-1); 
        dataPacket = fread(fid,21,'*double');
        
        INUdata.roll(frIDX+1)   = dataPacket(4); %degrees
        INUdata.pitch(frIDX+1)  = dataPacket(5); %degrees
        INUdata.azm(frIDX+1)    = dataPacket(6); %degrees
        INUdata.vroll(frIDX+1)  = (180/pi)*dataPacket(7); %convert to degrees/sec
        INUdata.vpitch(frIDX+1) = (180/pi)*dataPacket(8); %convert to degrees/sec
        INUdata.vazm(frIDX+1)   = (180/pi)*dataPacket(9); %convert to degrees/sec
        INUdata.aroll(frIDX+1)  = dataPacket(10); %meters/sec^2
        INUdata.apitch(frIDX+1) = dataPacket(11); %meters/sec^2
        INUdata.aazm(frIDX+1)   = dataPacket(12); %meters/sec^2
    elseif watson
        %seek past header and timestamps to get to 17-byte data packets
        fseek(fid,byteOffset+17*8*f,-1);
        dataPacket = fread(fid,17,'*double');
        
        INUdata.omega(frIDX+1) = dataPacket(14); %degrees
        INUdata.phi(frIDX+1)   = dataPacket(13); %degrees
        
        %subtract from 360deg to look like streams5 output
        INUdata.kappa(frIDX+1) = 360 - dataPacket(15); %degrees
        INUdata.vel(frIDX+1)   = dataPacket(16); %units? data values correct?
    end
end
fclose(fid);