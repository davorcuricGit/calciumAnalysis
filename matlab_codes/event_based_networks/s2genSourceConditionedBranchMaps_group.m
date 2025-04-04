clear all
setpathsloc = '/scratch/davorcuric/WFOIanalysisfunctions/';

run([setpathsloc 'setpaths.m'])
%dorsalMapLoc = [setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'];
load([setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'],'dorsalMaps');
load('allRecordings.mat')
expLoc = '/scratch/davorcuric/Sources_sinks_analysis/';

dmask = load('/scratch/davorcuric/maskDavor.mat'); dmask = dmask.Mask_Davor;
dmask = spatialBlockDownsample(double(dmask), 2, false); dmask(dmask ~= 0) = 1;
validPixels = find(dmask == 1);
conds = uniquecell(T.condition);
%% get the recordings we want

condsindexlist = [39];% 30 41];

for ci = condsindexlist;
    numRecs = [];
    numRecs = [find(strcmp(T.condition,conds{ci}))];


    warp = 1;
    err = 0;
    segmentDurationThreshold = 500;
    HKradius = 8;

    thresh = 1;

    clf

    downSample = 2;
    sz = [256 256]/downSample; clear downSample

    singleRoot = true;

    rootnum = 1; %get the single roots


  
    for r1 = [1:64]
        r2 = r1-1;
        clear profile_term rootProfile_term numAvs
        rootProfile_term = [];
        rootProfile_act = [];
        numAvs = [];
        count2 = 1;

        for iexp = numRecs'
                recLoc = T(iexp,:).paths{1};
    recName = T(iexp,:).names{1};
    recDataset = T(iexp,:).dataset{1};
    recCondition = T(iexp,:).condition{1};
    recMouse = T(iexp,:).mouse{1};
    recID = T(iexp,:).rec_id(1);


            [iexp]
                
            %skip outlier recordings    
            if T(iexp,:).outlier
                continue
            end

            %get the allen parcellation
            recLoc = T(iexp,:).paths{1};
            load([ recLoc '/transform.mat'])
            mask = output_warped(mask,tform,dorsalMaps);
            mask = spatialBlockDownsample(mask, 2, false);
            
            mask = mask.*dmask;

            allenMask = mask;



            RSCL = find(mask == r1);% | mask == 55);
            if isempty(RSCL);
                continue;
            end


            condmask = mask;
            condmask(condmask ~= r1) = 0;% & condmask ~= r2) = 0;
            condmask(condmask ~= 0) = 3;

            %get the avalanches
            
            loaddir = [expLoc '/results/HKAvalanches/HKradius= ' num2str(HKradius) '/imwarp/' ...
            recDataset '/' recCondition '/' 'recid=' num2str(recID) '/Thresh=' num2str(thresh) '/'];
            load([loaddir 'HKAvs.mat'])

            SS_source = resultToSave{1};
            DD = resultToSave{2};
            merged = resultToSave{3};
            source = resultToSave{4};
            sink = resultToSave{6};

            source = [source{:}];
            sink = [sink{:}];
            merged = [merged{:}];

            %take only the single roots
            idx = find(merged == 0);
            source = {source{idx}};
            sink = {sink{idx}};


            validPixels = resultToSave{7};
            mask = embeddIntoFOV(validPixels, validPixels, sz);

            S = [];
            count = 1;
            profile_term = zeros(1, length(validPixels));
            profile_act = zeros(1, length(validPixels));
            for i = 1:length(source)

                rr = source{i};


                %find the center of the source
                temp = zeros(size(validPixels));
                temp(rr) = 1;
                temp = embeddIntoFOV(temp, validPixels, [128 128]);
                stats = regionprops(temp);
                centroid = arrayfun(@round, stats.Centroid);
                centroid = sub2ind([128,128],centroid(2), centroid(1));



                if ~ismember(centroid, [RSCL])
                    continue
                end
                profile_term(count,sink{i}) = 1; %counts the number of terminating pixels
                profile_act(count, sink{i}) = 1/length(sink{i}); %counts the number of activations
                count = count + 1;
            end
 
  


            if ~exist('profile_term', 'var')
                stop
                continue
            end

            profile_term = embeddIntoFOV(profile_term', validPixels, sz);
            profile_act = embeddIntoFOV(profile_act', validPixels, sz);

            mask = embeddIntoFOV(validPixels, validPixels, sz);
            mask = abs(mask - 1);
            invertIdx = find(mask == 1);
            profile_term(invertIdx) = nan;
            profile_act(invertIdx) = nan;

            rootProfile_term(:,:,count2 ) = sum(profile_term, 3);
            rootProfile_act(:,:,count2 ) = sum(profile_act, 3);
            numAvs(:,:,count2 ) = size(profile_term,3);
            count2 = count2 + 1;
        end


        if isempty(rootProfile_term)
          continue
        end
        
        resultToSave = {rootProfile_term, rootProfile_act, numAvs};
        savedir = [expLoc '/results/conditionalsinks/matfiles/' conds{ci} '/thresh=' num2str(thresh) '/'];
        checkifDirExists(savedir)
        save([savedir 'sink_conditionalon=' num2str(r1) '_' labels{r1} '.mat'], 'resultToSave');
        %end
% 
%    
%         outline = allenMask;
%         outline(outline == round(outline)) = 0;
%         outline(outline ~= 0) = 1;
% 
%         condcenter = condmask;
%         condcenter(condcenter > 0) = 1;
%         stats = regionprops(condcenter);
%         try
%             condcenter = arrayfun(@round, stats.Centroid);
%         catch
%             condcenter = [0,0];
%         end
%         %condcenter = sub2ind([128,128],condcenter(2), condcenter(1));       
% 
%         imagesc(log10(sum(rootProfile_term,3) + 0.1*outline))
%         hold on
%         plot(condcenter(1),condcenter(2),'ro')
%         title(['conditioned on ' labels{r1}])
% 
%         savedir = ['/home/dcuric/Documents/' 'sourcesinksanalysis/results/conditionalsinks/figures/' conds{ci} '/'];
%         checkifDirExists(savedir)
%         saveas(gcf,[savedir 'sink_conditionalon=' num2str(r1) '_' labels{r1} '.png']);
    end
end
%%



