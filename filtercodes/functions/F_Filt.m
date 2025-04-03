
% Filters 2d image stack
% I is X,Y,t image stack
% Wn is filter window (Hz)
% Fs is sampling rate (Hz)
% parallel flag 1 to use parfor loop
% parallel flag 0 to use for loop

function ImgF=F_Filt(I,Wn,Fs,parallel)
% I: input 3D images
tic

if nargin <2; Wn = [0.5 6]; end %default bandpass filter 0.5 to 6 Hz
if nargin <3; Fs=150; end %default sample rate 150Hz
if nargin <4; parallel=0; end %default parallel off

[Sx,Sy,Sz]=size(I);

f = fspecial('disk',1);I= imfilter(I,f,'same');
n = 3; r = 0.5; ftype = 'bandpass';
[b,a]=cheby1(n,r,Wn/(Fs/2));%define filtering chebyl
I=reshape(I,Sx*Sy,Sz);ImgF=single(I);
if parallel
parfor i=1:Sx*Sy;ImgF(i,:)=single(filtfilt(b,a,double(I(i,:))));end % filtfilt
else;for i=1:Sx*Sy;ImgF(i,:)=single(filtfilt(b,a,double(I(i,:))));end % filtfilt
end
ImgF=reshape(ImgF,Sx,Sy,Sz);

disp(['Finshed Filt program in ' num2str(round(toc)) 's']);
end
