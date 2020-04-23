function [Time,info] = readVLTime(filepath)
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

%% Read Raw Header and Validate Requested Frame Numbers
fid=fopen(filepath);
if fid==-1, error(['Cannot open file: ' filepath]); end        
ncols               = fread(fid,1,'uint32');
nrows               = fread(fid,1,'uint32');
pixDepth            = fread(fid,1,'uint32');
nRawFrames          = fread(fid,1,'uint64');
byteOffsetToImagery = fread(fid,1,'uint64'); %nbytes of header + timestamps

VLFrNum = 1:nRawFrames;
nVLFr = numel(VLFrNum);

%% Read Raw Timestamps
%seek from beginning of file (bof) past 28byte header to 8byte timestamps
%to the first requested raw frame index
fseek(fid,28,'bof');

%read raw 8byte timestamps to double 
epochTimes  = fread(fid,nVLFr,'uint64=>double');

%convert from Epoch time (100nsec intervals since Jan 1, 1601) 
ymdhms = datevec(datenum(1601,1,1,0,0,0) + epochTimes/86400/1e7);
Time = ymdhms;
    
fclose(fid);

%% Output

    info.FilePath       = filepath;
    info.YMDHMS         = ymdhms;
    info.Time           = datetime(ymdhms);
    info.VLFrNum        = VLFrNum';
    info.VLRawHeader.nRows               = nrows;
    info.VLRawHeader.nCols               = ncols;
    info.VLRawHeader.PixDepth            = pixDepth;
    info.VLRawHeader.ByteOffsetToImagery = byteOffsetToImagery;

