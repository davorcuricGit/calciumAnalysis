


parent='Z:\acute_stress';
f=dir(fullfile(parent,'**','*fChan.dat'));

for n=1:length(f)
	
	folderdest=strrep(f(n).folder,'Z:\','H:\');
	
	mkdir(folderdest);
	disp([int2str(n) ': ' folderdest])

	copyfile(fullfile(f(n).folder,f(n).name),folderdest)
	copyfile(fullfile(f(n).folder,'transform.mat'),folderdest)
	copyfile(fullfile(f(n).folder,'Data_Fluo.mat'),folderdest)
	copyfile(fullfile(f(n).folder,'frame_error.mat'),folderdest)

end

