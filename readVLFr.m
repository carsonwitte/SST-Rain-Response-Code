function [fr,info] = readVLFr(filepath,VLFrNum,nFramesPerInc,rawFrameNumOffset,varargin)
% [fr info] = readVLFr(filepath,VLFrNum,nFramesPerInc,rawFrameNumOffset,varargin)
% 
% Reads an image or stack of images from a raw binary file.  This 
% raw file should be exported from a Streams VL Archive file and should 
% contain the standard 28-byte Streams imagery export header, the image 
% timestamps, and the image data at the end of the file.
%
% This is meant to replace readFr.m, with primary differences:
%    *returns images as cols x rows, where readFr returned rows x cols
%    *returns a 3D frame stack, instead of only 1 frame per call, with no
%     need to pre-initialized the stack variable/memory outside of the call
%    *checks RAM to insure there is enough available for requested stack
%    *returns 30% - 60% faster than readFr.m
%    *uses the stand-alone function plotVLFr.m to plot a single frame or 
%     a movie of a frame stack (see VARARGIN description below).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SYNTAX: fr = readVLFr(filepath,readFr)
%         fr = readVLFr(filepath,readFr,[],[],'PLOT',varargin)
%         fr = readVLFr(filepath,readFr,[],rawFrameNumOffset,'PLOT',varargin)
%         fr = readVLFr(filepath,readFr,nFramesPerInc,[],'PLOT',varargin)
%         fr = readVLFr(filepath,readFr,nFramesPerInc,rawFrameNumOffset,'PLOT',varargin)
%        [fr,info] = readVLFr( __ )
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs: FILEPATH 
%            A char vector, the full path location of the file.
% 
%         VLFRNUM   
%            A numeric vector of VL frame numbers to read from the RAW 
%            file.  Frame numbers must be increasing and equally
%            spaced/incremented.  
%
%            For example, VLFRNUM = 0:3:10 will return the four VL frames 
%            0,3,6,9.  NOTE in this example the increment between frames is
%            three, and the optional NFRAMESPERINC can be used to return
%            more than one frame per frame increment (see below).
% 
%         NFRAMESPERINC (optional)
%            A single number option always greater than 1 and less than 
%            the equally spaced frame increments of VLFRNUM.  This value 
%            allows the return of more than 1 image per frame increment.
% 
%            For example, VLFRNUM = 0:3:10 w/ option NFRAMESPERINC = 2 
%            will return the eight VL frames 0,1,3,4,6,7,9,10.
% 
%         RAWFRAMENUMOFFSET (when raw file does not contain all VL frames)
%            A single number to denote the first VL frame number of 
%            the current raw file.  A value should be entered if the 
%            first frame of the raw file is NOT frame 0 of the archive VL.
%            Otherwise VLFRNUM will be referenced incorrectly and not
%            agree with frame numbers in the VL Archive file.  
%
%            NOTE the only reason to export raw files that are only a 
%            portion of the VL Archive files is to save space.  The raw 
%            offset value is unique for each device in a given VL.  The 
%            value should be denoted in the raw filename when the raw file
%            is exported.  This is useful in cases where you want to test
%            image processing locally on a small number of exported frames, 
%            or in cases where VL devices are recording meaningless data
%            (such as when a UAV payload is readying for launch).
%
%         VARARGIN
%            When 'PLOT' is entered as the first arg, plotVLFr.m is called.
%            plotVLFr.m has its own arg options to handle single images or 
%            a stack of images.  See plotVLFr.m for all available options.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs: FR   
%            UINT16 image data of size MxNxnFrames where MxN is Rows x Cols
%            and nFrames is the total number of requested frames.
% 
%          INFO 
%            A struct w/ frame times (usually UTC unless Streams was not
%            sync'd by GPS), frame numbers, raw offset (if any), and the
%            raw export header information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Scott Brown
% 9 Jul 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input Arg Control
if nargin<2
    error('You must enter at least the first two arguments.')
elseif nargin==2
    rawFrameNumOffset = 0;
    nFramesPerInc     = 1;
elseif nargin==3
    rawFrameNumOffset = 0;
end

if isempty(rawFrameNumOffset), rawFrameNumOffset=0; end
if isempty(nFramesPerInc), nFramesPerInc=1; end

%% Determine & Validate Frame Increment
frdiff = diff(VLFrNum);   
if isempty(frdiff) || numel(VLFrNum)==1
    frameInc = 1;
elseif ~all(frdiff==frdiff(1)) || frdiff(1)<1
    error(['VLFRNUM should be an equally spaced, increasing ',...
           'numeric vector.'])
elseif numel(VLFrNum)>=2    
    frameInc = frdiff(1);
end

%% Add Additional Requested VL Frame Numbers if NFRAMESPERINC > 1
nVLFr = numel(VLFrNum);
if nFramesPerInc>1 && nFramesPerInc>=frameInc
    error(['Cannot return nframes-per-increment greater than or equal ',...
           'to the frame increment value.  Decrease the requested ',...
           'NFRAMESPERINC, or increase the frame increment.'])
elseif nFramesPerInc>1    
    keepFrNumyn = repmat([true(1,nFramesPerInc) false(1,frameInc-nFramesPerInc-1)],...
                         [1 nVLFr]);

    tmpFrNum = VLFrNum(1):VLFrNum(end);
    addlIDX  = repmat(0:nVLFr-1,[nFramesPerInc 1]);
    addlIDX  = addlIDX(:)';
    VLFrNum  = tmpFrNum(keepFrNumyn)+addlIDX;
end
%update number of read frames
nVLFr = numel(VLFrNum);

%% Read Raw Header and Validate Requested Frame Numbers
fid=fopen(filepath);
if fid==-1, error(['Cannot open file: ' filepath]); end        
ncols               = fread(fid,1,'uint32');
nrows               = fread(fid,1,'uint32');
pixDepth            = fread(fid,1,'uint32');
nRawFrames          = fread(fid,1,'uint64');
byteOffsetToImagery = fread(fid,1,'uint64'); %nbytes of header + timestamps

%check requested VL frame numbers to make sure they exist in this RAW file
nTotalVLFrames = nRawFrames+rawFrameNumOffset;
if any(VLFrNum > nTotalVLFrames-1) || any(VLFrNum < rawFrameNumOffset)   
    fclose(fid);
    error(['VL frame numbers available in this RAW file are ',...
           num2str(rawFrameNumOffset) '-',... 
           num2str(nRawFrames+rawFrameNumOffset-1) '.'])
end
nFramePixels = ncols*nrows;

%% Warning if Requesting 75-99% of Memory (RAM), Error if > 99%
if ispc()
    [~,sys]     = memory();
    unusedRAMGB = sys.PhysicalMemory.Available/1000^3;
elseif ismac()
    unusedRAMGB = system('sysctl -n hw.physmem | awk ''{print $2}''');
elseif isunix()
    [~,w] = unix('free | grep Mem');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    unusedRAMGB = stats(3);
end
%2bytes per pixel x nFrames convert bytes to GB
sizeAllFramesGB = nFramePixels*nVLFr*2/1000^3;
usedPrct        = sizeAllFramesGB/unusedRAMGB*100;

if usedPrct>75 && usedPrct<100
    warning(['Number of requested frames at size ' num2str(ncols) 'x',...
              num2str(nrows) ' requires ' num2str(sizeAllFramesGB),...
             'GB of physical memory (RAM), which is currently ',...
              sprintf('%2.1f',usedPrct) '% of total available.'])
elseif usedPrct>=100
    fclose(fid);
    error(['Number of requested frames at size ' num2str(ncols) 'x',...
            num2str(nrows) ' requires ' num2str(sizeAllFramesGB),...
           'GB, but there is only ' num2str(unusedRAMGB) 'GB of total ',...
           'available physical memory (RAM)'])
end

%% Read Raw Timestamps
%seek from beginning of file (bof) past 28byte header to 8byte timestamps
%to the first requested raw frame index
sRawByteIDX = VLFrNum(1)-rawFrameNumOffset;
fseek(fid,28+8*sRawByteIDX,'bof');
%for frameInc>1, skip this many bytes between frames
skip8byte = 8*(frameInc-nFramesPerInc); 
%read raw 8byte timestamps to double 
epochTimes  = fread(fid,nVLFr,...
                       [num2str(nFramesPerInc) '*uint64=>double'],...
                       skip8byte);
%convert from Epoch time (100nsec intervals since Jan 1, 1601) 
ymdhms = datevec(datenum(1601,1,1,0,0,0) + epochTimes/86400/1e7);

%% Read Raw Imagery Data
%seek from beginning of file (bof) past 28byte header 8byte timestamps to
%the first requested raw frame index
fseek(fid,byteOffsetToImagery + 2*nFramePixels*sRawByteIDX,'bof');
%for each frame increment, move 2bytes*image*[diff between frame increments]
skip2byte = 2*nFramePixels*(frameInc-nFramesPerInc); 
%read raw 2byte pixel data to uint16
frvector  = fread(fid,nFramePixels*nVLFr,...
                     [num2str(nFramePixels*nFramesPerInc) '*uint16=>uint16'],...
                     skip2byte);
%reshape to a 3D stack of cols x rows x nFrames
fr = reshape(frvector,[ncols nrows nVLFr]);
    
fclose(fid);

%% Return Info on Requested Frame(s), also VL Metadata
if nargout==2
    info.FilePath       = filepath;
    info.YMDHMS         = ymdhms;
    info.Time           = datetime(ymdhms);
    info.VLFrNum        = VLFrNum';
    info.RawFrNumOffset = rawFrameNumOffset;
    info.nRawFrames     = nRawFrames;
    info.nTotalVLFrames = nTotalVLFrames;
    info.VLRawHeader.nRows               = nrows;
    info.VLRawHeader.nCols               = ncols;
    info.VLRawHeader.PixDepth            = pixDepth;
    info.VLRawHeader.ByteOffsetToImagery = byteOffsetToImagery;
end

%% Plotting Will Occur Only if First VARARGIN Arg is 'PLOT'
if numel(varargin)>0 && strcmpi(varargin{1},'PLOT')
    plotVLFr(fr,info,[],varargin{2:end})
end


