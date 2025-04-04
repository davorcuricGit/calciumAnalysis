function F_SaveRAW(Data,filename)
tic
[Sx,Sy,Sz]=size(Data);
if Sz<1; error(['F_SaveRAW;The Data input is 2D ' num2str(Sx) ' ' num2str(Sy) ' ' num2str(Sz) '.']);end
for i=1:Sz;Data(:,:,i)=Data(:,:,i)';end %flip x y
filename=[filename '_float32_' num2str(Sx) '-' num2str(Sy) '-' num2str(Sz) '.raw'];
%disp(filename);
filename
fid = fopen(filename,'w', 'b');
fid
fwrite(fid, Data, 'float32');
fclose(fid);

disp(['Write file finised in ' num2str(round(toc)) 's: ' filename]);