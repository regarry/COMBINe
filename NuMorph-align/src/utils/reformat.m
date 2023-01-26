function translation_table = reformat(config)
% save translation parameters from alignment_table --> for python loadinng

load(fullfile(config.output_directory,'variables','alignment_table.mat'),...
    'alignment_table');
t_size = size(alignment_table);
for i = 2:length(config.markers)
    if i>1 && any(~ismember(config.align_channels,i))
        continue
    end
    
    c = char(config.markers(i));
    params = cell(t_size);
    for j = 1:t_size(1)
        for k = 1:t_size(2)
            X_disp = alignment_table{j,k}.(strcat('X_Shift_',c));
            Y_disp = alignment_table{j,k}.(strcat('Y_Shift_',c));
            params{j,k} = cat(2, X_disp, Y_disp);
        end
    end
    translation_table.(config.markers(i)) = params;
    save_path = fullfile(config.output_directory,'variables','translation_table.mat');
    save(save_path,'translation_table')
    
end

