%this code requires that you have run branchmaps first
clear all
setpathsloc = '/scratch/davorcuric/WFOIanalysisfunctions/';

run([setpathsloc 'setpaths.m'])
%dorsalMapLoc = [setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'];
load([setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'],'dorsalMaps');
load('/scratch/davorcuric/allRecordings.mat')
expLoc = '/scratch/davorcuric/Sources_sinks_analysis/';

dmask = load('/scratch/davorcuric/maskDavor.mat'); dmask = dmask.Mask_Davor;
dmask = spatialBlockDownsample(double(dmask), 2, false); dmask(dmask ~= 0) = 1;
validPixels = find(dmask == 1);

load('/scratch/davorcuric/centroids.mat')
numpix = centroid{2};
nodepos = centroid{1};

% get the recordings we want
conds = uniquecell(T.condition);
condsindexlist = 41%[24 27 26];
numRecs = [];

numRecs = find(strcmp(T.condition,conds{condsindexlist(1)}));



%load this just to get the labels
recLoc = T(1,:).paths{1};
load([ recLoc '/transform.mat'])
clear mask

%
%numRecs = %height(T);

warp = 1;
err = 0;
segmentDurationThreshold = 500;
HKradius = 8;

thresh = 1;

%clf

downSample = 2;
sz = [256 256]/downSample;

xmin = 0.0033;

singleRoot = true;

rootnum = 1; %get the single roots

Adj_term = zeros(64,64);
Adj_act = zeros(64,64);
numActivations = cell(64,1);
totnumavs = [];
for r1 = [1:64]
    
    numActivations{r1} = 0;
    npix{r1} = numpix(r1); %need to recast as cell because of matlab table 6__6

        spatialmapname = [expLoc '/results/conditionalsinks/matfiles/' conds{condsindexlist(1)} '/' 'thresh=' num2str(thresh) '/'];
        spatialmapname = [spatialmapname 'sink_conditionalon=' num2str(r1) '_' labels{r1} '.mat'];
        
        if exist(spatialmapname)
            load(spatialmapname)
        else
            continue
        end
        
        if numpix(r1) < 50;
            continue
        end

        spatialmap_term = resultToSave{1};
        spatialmap_act = resultToSave{2};
        numavs = squeeze(resultToSave{3});
        
        if isempty(spatialmap_term);
            continue
        end


    for iexp = 1:size(spatialmap_term,3);
      
        %get the network
        recLoc = T(numRecs(iexp),:).paths{1};
        load([ recLoc '/transform.mat'])
        mask = output_warped(mask,tform,dorsalMaps);
        mask = spatialBlockDownsample(mask, 2, false);
        
        for i = 1:1:max(max(mask));
            temp = mask.*dmask;
            
            temp(temp ~= i) = 0;
            temp(temp == i) = 1;

            %if the region is too small take it out of the possible targets
            if sum(sum(temp)) < 50;
                continue
            end
            

            %Weights are proportional to the number of termination pixels
            Adj_term(r1, i) = Adj_term(r1, i) + (nansum(nansum(spatialmap_term(:,:,iexp).*temp)));
            
            %Weights are proportional to the number of avalanches
            Adj_act(r1, i) = Adj_act(r1, i) + (nansum(nansum(spatialmap_act(:,:,iexp).*temp)));


        end
        numActivations{r1} = numActivations{r1} + numavs(iexp);
        
    end
end

G_term = digraph(Adj_term/length(numRecs), labels);
G_act = digraph(Adj_act/length(numRecs), labels);

%

savedir = ['/home/dcuric/Documents/sourcesinksanalysis' '/results/graphs/' conds{condsindexlist(1)} '/' ];
checkifDirExists(savedir)

clear nodes
nodes = table();
nodes(:,'names') = labels;
nodes(:,'posx') =  num2cell(nodepos(:,1));
nodes(:,'posy') =  num2cell(nodepos(:,2));
nodes(:, 'area') = npix';
nodes(:, 'activations') = numActivations;
writetable(nodes, [savedir 'nodes.csv'])

edges_term = table();
edges_term(:,'source') = G_term.Edges.EndNodes(:,1);
edges_term(:,'target') = G_term.Edges.EndNodes(:,2);
edges_term(:, 'weight') = num2cell(G_term.Edges.Weight);
writetable(edges_term, [savedir 'edges_numTermination_thresh=' num2str(thresh) '.csv'])

edges_act = table();
edges_act(:,'source') = G_act.Edges.EndNodes(:,1);
edges_act(:,'target') = G_act.Edges.EndNodes(:,2);
edges_act(:, 'weight') = num2cell(G_act.Edges.Weight);
writetable(edges_act, [savedir 'edges_numActivation_thresh=' num2str(thresh) '.csv'])


% 'saved'
% 
% savedir2 = ['/home/dcuric/Documents/sourcesinksanalysis/results/graphs/' conds{condsindexlist(1)} '/' ];
% checkifDirExists(savedir2)
% save([savedir2 'graph.mat'], 'G')
%%
clf
plot(G_term.Edges.Weight./G_act.Edges.Weight)
%%
clf
plot(sum(Adj_act))
hold on
plot(cell2mat(numActivations))