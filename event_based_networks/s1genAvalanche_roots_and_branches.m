
clear all
setpathsloc = '/scratch/davorcuric/WFOIanalysisfunctions/';

run([setpathsloc 'setpaths.m'])
dorsalMapLoc = [setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'];
load([setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'],'dorsalMaps');
clear dorsalMapLoc
load('allRecordings.mat')
load([setpathsloc '/functions/graphTheory/pariwiseDistanceMatrix_sz=128x128.mat'])
expLoc = '/scratch/davorcuric/Sources_sinks_analysis/';
conds = uniquecell(T.condition);
%%

%reclist = find(contains(T.condition,'ket100') | contains(T.dataset,'GCamp6_hip_rosbengal'));
reclist = find(contains(T.dataset, 'seizure'));

conditionslist = {'Awake', 'Iso1', 'Iso2', 'Ket10', 'ket100', 'pento12.5', 'Pent80_30'};
reclist = [];
for i = 1:length(conditionslist)
    reclist = [reclist find(contains(T.condition, conditionslist{i}))']
end


%% because we use a universal mask we can get the network right away since the distance matrix takes up a lot of memory

HKradius = 8;

    %get the mask
    mask = getfield(load('maskDavor.mat', 'Mask_Davor'), 'Mask_Davor');
    mask = spatialBlockDownsample(single(mask/255), 2);
    mask(mask ~= 1) = NaN;
    validPixels = single(find(mask == 1));
    
    'gen network'
    adjmat = D(validPixels, validPixels);
    adjmat(adjmat > HKradius) = 0;
    adjmat(adjmat ~= 0) = 1;
    clear network
    for i = 1:size(adjmat,1)
        network{i} = single(find(adjmat(i,:) == 1));
        network{i} = [network{i} i];
    end

    clear D

%%
clear ImgF

numRecs = height(T);

warp = 1;
err = 0;



segmentDurationThreshold = 500;

threshvals = [1:3];

%The image warping can cause a lot of memory use. batchblocks is the number
%of blocks the recording is split into, warped, and recombined into.
batchblocks = 4;

for iexp = reclist(2:end)
    recLoc = T(iexp,:).paths{1};
    recName = T(iexp,:).names{1};
    recDataset = T(iexp,:).dataset{1};
    recCondition = T(iexp,:).condition{1};
    recMouse = T(iexp,:).mouse{1};
    recID = T(iexp,:).rec_id(1);
    

    tSteps = T(iexp,:).durations{1};%project.recordings{iexp}.recordingDuration;
    [segmentToKeep,  goodFrames] = getSegmentToKeep(T(iexp,:).frameoffset{1}, tSteps, [recLoc '/' recName '.raw'], batchblocks, segmentDurationThreshold );

    if isempty(segmentToKeep); continue; end

   % 
    %load the data
    ImgF = F_ReadRAW([recLoc '/' recName], [256, 256, tSteps], 'float32', warp, err, 0, batchblocks);
 %


    clear tform tform_old points coordinates motionIndex

    'downsampling'
    tic()
    ImgF = spatialBlockDownsample(ImgF, 2);
    ImgF = zscore_independent(ImgF);
    
    ImgF = ImgF.*mask;
%
    for ii = 1:length(threshvals)%4:6]%[1:.25:3]
        thresh = threshvals(ii);
        count = 1;
%
        clear SS DD labeledFrame roots rootTimes branches branchTimes merged S D
        for i = 1:length(segmentToKeep)
            segTime = goodFrames{segmentToKeep(i)};
            seg = ImgF(:,:,segTime);

            seg(seg < thresh) = 0;
            seg(seg~= 0) = 1;

            if sum(sum(seg)) == 0; continue; end
            seg = reshape(seg, size(seg,1)*size(seg,2), []);
            seg = seg(validPixels, :);

            %get rid of any zero activity periods
            trace = (sum(seg));
            trace(1) = 0; trace(end) = 0;
            if sum(trace) == 0; continue; end
            [~, segTimes] = segmentSeries(trace);
            segTimes = {segTimes{2:2:end}};

            tic()

            for j = 1:length(segTimes)
                [iexp i/length(segmentToKeep) j j/length(segTimes) ]
                T0 = segTimes{j}(1);
                [S{count} D{count}, merged{count}, ~, roots{count}, rootTimes{count}, branches{count}] = getAvalanches(seg(:, segTimes{j}), network, adjmat, validPixels, T0);
                
                %%flip reverses the segment so branches become roots
                %[~, ~, ~, branches{count}, branchTimes{count}] = getAvalanches(flip(seg(:, segTimes{j}),2), network, adjmat, validPixels, T0);
                count = count + 1;

                
            end
            toc()
        end
        if exist('S') ~= 1;
            break;
        end
        
        savedir = [expLoc '/results/HKAvalanches/HKradius= ' num2str(HKradius) '/imwarp/' ...
            recDataset '/' recCondition '/' 'recid=' num2str(recID) '/Thresh=' num2str(thresh) '/'];

        savename = 'HKAvs';
        checkifDirExists(savedir)
        resultToSave = {S, D, merged, roots, rootTimes, branches, validPixels};
        clear SS DD roots rootTimes branchTimes branchTimes merged branches S D

        save([savedir savename '.mat'], 'resultToSave')

    end

    clear ImgF ImgF2
end






function [segmentToKeep, goodFrames] = getSegmentToKeep(T0, recLen, recordingName, batchblocks, segmentDurationThreshold );
%      get the frames with no excessive motion
%     some of the recordings are split up. This can cause problems with the motionIndex file.
%     the stupid way around this is to force an error which will then return the length of the motionIndex file,
%     which is then seperated



tSteps = 1000000;
try 
motionIndex = fetchMotionErrorFile( recordingName, [256, 256, tSteps]);
catch exception
    if strcmp(exception.identifier, 'MATLAB:badsubscript')
        tSteps = exception.message;
       tSteps = strsplit(tSteps, ' ');
       tSteps = str2num(tSteps{end}(1:end-1));
       motionIndex = fetchMotionErrorFile( recordingName, [256, 256, tSteps]);
    end
end

    %find segments of recording we want to keep (ie. long with no motion
    %artifact. If no segments are found then continue to next animal
    accTime = 1:tSteps; 
    
    accTime = accTime(T0:T0+recLen-1).*motionIndex(T0:T0+recLen-1)';
    accTime = accTime(1:end - mod(recLen, batchblocks)); %batchblocks may remove a few frames from the end
    
    %now make sure we are counting from one
    accTime = accTime - T0 + 1;
    accTime(accTime < 0) = 0;

    s = segmentSeries(accTime);
    goodFrames = s(2:2:end);
    %goodFrames = [goodFrames{:}];
    badFrames = s(1:2:end);
    Lgood = cellfun(@length, goodFrames);
    Lbad = cellfun(@length, badFrames);
    segmentToKeep = find(Lgood > segmentDurationThreshold);
end
