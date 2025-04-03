% Read .dat file into a good format
% requires a folder containing a single Data_Fluo.mat file
% and a single fChan.dat file

%returns the image stack in an x,y,z array
%also returns the sampling rate Fs

% April 2021 edit: flag for different channels:
% 0 default/null: one channel
% 1 567 fluo
% 2 475 fluo
% 3 red
% 4 green

function [I,Fs]=F_ReadDAT(folder,maxFrames,CHAN_FLAG)
% Load fChan file
if nargin<3
 CHAN_FLAG=0;
end

switch CHAN_FLAG
 case 0
  disp('default chosen..')
  mat_file='Data_Fluo.mat';
  dat_file='fChan.dat';
 case 1
  disp('Fluo 567 chosen..')
  mat_file='Data_Fluo_567.mat';
  dat_file='fChan_567.dat';
 case 2
  disp('Fluo 475 chosen..')
  mat_file='Data_Fluo_475.mat';
  dat_file='fChan_475.dat';
 case 3
  disp('Red chosen..')
  mat_file='Data_red.mat';
  dat_file='rChan.dat';
 case 4
  disp('Green chosen..')
  mat_file='Data_green.mat';
  dat_file='gChan.dat';
 otherwise
  error('not a valid channel selection..')  
end

tic
Info = matfile(fullfile(folder,mat_file));
Sy=Info.datSize(1,1);
Sx=Info.datSize(1,2);
Sz=Info.datLength; %number of frames
Fs=Info.Freq; %sampling rate

if nargin>1 %if more than 1 input
if maxFrames<Sz&&maxFrames>0 %if nFrames is greater than maxFrames
Sz=maxFrames; %read only maxFrames
end
end

if Fs==0
Fs=30;
end

filename=fullfile(folder,dat_file);
I=repmat(single(0),[Sx,Sy,Sz]);
fid = fopen(filename);
for m=1:Sz
dat = fread(fid,[Sx,Sy],'single');
I(:,:,m)=dat';
end
fclose(fid);

disp(['Read file finised in ' num2str(round(toc)) 's: ' filename]);