clear all
run('/home/dcuric/Documents/calciumAnalysis/matlab_codes/init_analysis.m')

load(['/home/dcuric/Documents/calciumAnalysis/matlab_codes/allRecordings_carnot.mat']);
T = project.project_table;
project_root = project.project_root;
conds = uniquecell(T.condition);

%%

reclist = T;%(find(contains(T.condition,'Awake')),:);% | (contains(T.condition,'Iso1')));

%reclist = find(contains(T.dataset,'GCamp6_hip_rosbengal') | contains(T.dataset,'seizure') | contains(T.condition,'ket100'));

%reclist = find(contains(T.dataset,'accute_stress'));


%%
[validPixels,dmask] = get_universal_mask(setpathsloc);
    

%%

clear ImgF

%The image warping can cause a lot of memory use. batchblocks is the number
%of blocks the recording is split into, warped, and recombined into.
batchblocks = 4;

FCderivativeTable = table(); 


for iexp = 1:reclist'
    row = reclist(iexp,:);

    recording_path = [project_root '/raw_data/' row.paths{1}];
    recording_name = [recording_path '/' row.names{1}];

    [ImgF, good_frames,segment_to_keep] = get_calcium_recording(recording_name,row, 'down_sample', 1, 'batch_blocks', batchblocks);
    ImgF = ImgF(:,:,[good_frames{segment_to_keep}]);
    ImgF = threeD_to_twoD(ImgF, 'validPixels', validPixels);

    
    subnetworks = get_subnetworks(recording_path, dorsalMaps, 'global_mask', dmask);
 

    params.len_control = 0;
    params.random_sample = 1;
    FCrow = make_FC(ImgF, subnetworks, row, project_root, params);
    FCderivativeTable = [FCderivativeTable; FCrow];
    
    %control the length
    L = 4000;
    if size(ImgF, 2) > L
    params.len_control = L;
    params.random_sample = 0;
    FCrow = make_FC(ImgF, subnetworks, row, project_root, params);
    FCderivativeTable = [FCderivativeTable; FCrow];
    
    params.len_control = L;
    params.random_sample = 1;
    FCrow = make_FC(ImgF, subnetworks, row, project_root, params);
    FCderivativeTable = [FCderivativeTable; FCrow];
    end
    
    save('FCderivativeTable.mat', 'FCderivativeTable');
    
end

