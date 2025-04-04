
function FCrow = make_FC(ImgF, subnetworks, row, project_root, params)
    %generate an FC from ImgF
    %ImgF is N x Time
    %subnetworks is cells of pixels corresponding to subnetworks
    %row is from the project table
    %outputs the derivative information adn saves FC

  clear v
    for i = 1:length(subnetworks);
        px = subnetworks{i};

        if params.random_sample;
            px = datasample(px, floor(length(px)/2), 'Replace', false);
        end

        if params.len_control > 0;
            v(i,:) = zscore(nanmean(ImgF(px,1:params.len_control)));
        else
            v(i,:) = zscore(nanmean(ImgF(px,:)));
        end
    end
    
    C = corrcoef(v');
    
    main_save_dir = ['/derivatives/' row.paths{1}];
    
    save_dir = [main_save_dir];
    if params.len_control >0
        save_dir = [save_dir '/len_control=' num2str(params.len_control) '/'];
    end
    if params.random_sample
        save_dir = [save_dir '/random_sample/'];
    end
 

    checkifDirExists([project_root save_dir])
    save_name = ['FC_' strrep(row.names{1}, '.raw', '.mat')];

    FCrow = removevars(row, {'names', 'paths', 'durations', 'frameoffset', 'machine_p', 'height', 'width'});
    FCrow.names = convertCharsToStrings(save_name);
    FCrow.paths = convertCharsToStrings(save_dir);
    
    FCrow.type = 'FC';
    FCrow.len_control = params.len_control;
    FCrow.random_sampled = params.random_sample;
    FCrow.OG_length = size(v,2);
    FCrow.numROI = size(v,1);
    

    save([project_root save_dir save_name], 'C')


end
