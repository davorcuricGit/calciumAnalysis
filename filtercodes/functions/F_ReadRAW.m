% revised 10/13/20 DA
% should be cross compatible with previous F_ReadRAW
% 
% Folder must contain a 'transform.mat' file for warping
% Folder must contain a 'frame_error.mat' file for error removal

% warp: 1 to output warped image stack
% rm_error: 1 to remove movement error frames
% nFrames: number of random frames to read into dataset (DA edit Nov23/21)
% memsave: image warping can take a lot of memory. If memsave > 1 then split the
% array into memsave number of blocks, warp, and recombine. If memsave <= 1 (default = 0), warp whole recording.

function [ImgF,frame_idx]=F_ReadRAW(filename,Sxyz,Ftype,warp,rm_error,nFrames, memsave)


tic
if nargin < 2 | Sxyz==0 
Sxyz=filename;
Sxyz(Sxyz<'0'|Sxyz>'9')=' ';
Sxyz=str2num(Sxyz);
Sxyz=Sxyz(end-2:end); %last three numbers need to be x,y,t
end  
if nargin < 3 | Ftype==0; Ftype = 'float32'; end  
if nargin < 4;warp=0;end
if nargin < 5;rm_error=0;end
%
frame_idx=1:Sxyz(3);

if rm_error==1
	[pathname,~]=fileparts(filename);
load([pathname '/frame_error.mat'],'error');
index=remove_errors(error);
%ImgF=ImgF(:,:,index);

frame_idx=frame_idx(index);
end

if nargin < 6 | nFrames==0
	frame_idx=frame_idx;
else
	index_rand=randi(length(frame_idx),nFrames,1);
	frame_idx=frame_idx(index_rand);
end

if nargin < 7
    memsave = 0;
end


ImgF=repmat(single(0),[Sxyz(1),Sxyz(2),length(frame_idx)]);
fid = fopen(filename,'r','b');
% for i=1:Sxyz(3)
%     Data=fread(fid,[Sxyz(1),Sxyz(2)],['*' Ftype]);
%     ImgF(:,:,i)=Data';
% end
% fclose(fid);
for i=1:length(frame_idx)
	ptr=(frame_idx(i)-1)*Sxyz(1)*Sxyz(2)*4;
	fseek(fid,ptr,'bof');
	Data=fread(fid,[Sxyz(1),Sxyz(2)],['*' Ftype]);
	ImgF(:,:,i)=Data';
end
fclose(fid);

if warp==1
[pathname,~]=fileparts(filename);
load([pathname '/transform.mat'],'tform');
load('/scratch/davorcuric/allenDorsalMap_donovan.mat','dorsalMaps');

if memsave <= 1
    ImgF=output_warped(ImgF,tform,dorsalMaps);
else
    nblocks = memsave;
    blocks = 1:Sxyz(3);
    blocks = reshape(blocks(1:end - mod(Sxyz(3), nblocks)), [], nblocks)';
    imgFwarped = [];
    for b = 1:nblocks
           ImgF2=output_warped(ImgF(:,:,blocks(b,:)),tform,dorsalMaps);
           imgFwarped = cat(3, imgFwarped, ImgF2);
    end
    clear ImgF
    ImgF = imgFwarped;
    clear imgFwarped
    %ImgF2=output_warped(ImgF(:,:,1:floor(size(ImgF,3)/2)),tform,dorsalMaps);
    %ImgF=output_warped(ImgF(:,:,floor(size(ImgF,3)/2)+1:end),tform,dorsalMaps);
    %ImgF = cat(3, ImgF2, ImgF);
end
end

disp(['Read file finised in ' num2str(round(toc)) 's: ' filename]);
end


%%
%warp
function img=output_warped(ImgF,tform,dorsalMaps)
% use tform to output transformed
ImgF_warped=imwarp(ImgF,tform,'OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));

%clip into 256*256 window (HARD CODED)
img=ImgF_warped;
yind=(151-127):(152+127);
xind=15:(15+255);
img=img(xind,yind,:);

end

%remove errors
function [index]=remove_errors(error)

ol=isoutlier(error);
ol(1:5)=0; %ensure first value is not an outlier
ol(end)=0; %ensure last value is not an outlier
cross=diff(ol);

x1=find(cross==1); %start of outlier
x2=find(cross==-1); %end of outlier
x1=x1-5; %5 frame buffer
x2=x2+5; %5 frame buffer

%shift to find stable segments
x2(2:end+1)=x2;
x2(1)=1;
%add end
x1(end+1)=length(ol);

dur=x1-x2;
ii=find(dur>30); %choose only segments longer than 1 sec 
ind=[];
for n=1:length(ii)
ind=cat(2,ind,x2(ii(n)):x1(ii(n)));
end

index=false(length(ol),1);
index(ind)=1;
end

