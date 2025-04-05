function [ImgF, good_frames,segment_to_keep] = get_calcium_recording(recording_path,row, varargin);

i_p = inputParser;
i_p.addRequired('recording_path', @ischar)
i_p.addRequired('row', @istable);
i_p.addOptional('warp',1);
i_p.addOptional('err', 0);
i_p.addOptional('batch_blocks', 1);
i_p.addOptional('size', [256,256]);
i_p.addOptional('machine_p', 'float32');
i_p.addOptional('down_sample', 2);
i_p.addOptional('zscore_flag', logical(1));
i_p.addOptional('segmentDurationThreshold', 500);
%i_p.addOptional('mask', nan);

i_p.parse(recording_path, row, varargin{:});



tSteps = row.durations{1};
offset = row.frameoffset{1};

% if isnan(i_p.Results.mask)
%     mask = ones(size);
% else
%     mask = i_p.Results.mask;
% end

%get the motion files. If there is too mcuh motion then we don't even
%bother retreievinge file
[segment_to_keep,  good_frames] = get_segments_to_keep(offset, tSteps, recording_path, i_p.Results.batch_blocks, i_p.Results.segmentDurationThreshold );
if isempty(segment_to_keep);
    ImgF = 'too much motion';
else

    %load the data
    ImgF = F_ReadRAW(recording_path, [i_p.Results.size(1), i_p.Results.size(2), tSteps], i_p.Results.machine_p, i_p.Results.warp, i_p.Results.err, 0, i_p.Results.batch_blocks);

    clear tform tform_old points coordinates motionIndex

    
    if i_p.Results.down_sample > 1
        'downsampling'
        ImgF = spatialBlockDownsample(ImgF, i_p.Results.down_sample);
%        mask = spatialBlockDownsample(mask, i_p.Results.down_sample, false);
    end
    if i_p.Results.zscore_flag
        ImgF = zscore_independent(ImgF);
    end


end