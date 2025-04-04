clear all
%% get the files

d = dir(fullfile('/scratch4/acute_stress/','**','*fChan.dat'));
folders = {d.folder};

band = [0.1 15];

for i = 1:length(folders)
    filedir = folders{i};
    datflou = load([filedir '/Data_Fluo.mat']);
    fStep = datflou.datLength;
    freq = datflou.Freq;
    

    [ImgF,Fs] = F_ReadDAT(filedir, fStep);

    F0 = mean(ImgF, 3);
    ImgF=(ImgF-F0)./F0*100;
    ImgF = F_Filt(ImgF, band, Fs);

    savename = [filedir '/Spontaneous_' num2str(freq) 'Hz_filt' ...
        strrep(num2str(band(1)), '.', '') '-' strrep(num2str(band(2)), '.', '')
        ];

    F_SaveRAW(ImgF, savename)
    delete([filedir '/fChan.dat'])
    clear ImgF
end

%%


%%


%%


F_SaveRAW(filtF, savename)
