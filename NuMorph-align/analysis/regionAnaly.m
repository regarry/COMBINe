% Region-wise quantitative analysis between groups

%% Regional quantitation 

dataset_dir = 'path\to\all\datasets\and\quantifications';
classes=["yellow_neuron","yellow_astrocyte","green_neuron",...
    "green_astrocyte","red_neuron","red_astrocyte"];
header={'x','y','z','xt','yt','zt','id','name','sub1','sub2','sub3'};
save_table = true;

result_dir = fullfile(dataset_dir, 'analysis');
mkdir(result_dir);
[data_subdirs,~] = GetSubDirsFirstLevelOnly(dataset_dir);

volumes = cell(1,length(data_subdirs));
counts = cell(1,length(data_subdirs));
densities = cell(1,length(data_subdirs));
df_results = cell(1,length(data_subdirs));
for i = 1:length(data_subdirs)
    data_dir = data_subdirs{i};
    [volumes{1,i},counts{1,i},densities{1,i},...
        df_results{1,i}] = quantify(data_dir,classes, header, save_table);
end

save(fullfile(result_dir,'results.mat'), 'volumes','counts','densities','df_results');

%% Assign groups and get statistics

groups = {[1,4,6],[2,3,5]}; % assign group indexes and divide into two groups
group_names = ["+/+","F/+"]; % names of groups

df_template = readtable(fullfile(fileparts(which('NM_config')),'annotations','structure_template.csv'));
indexes = df_template.index;

% seperate results by group
group_results = cell(1,length(groups));
for j = 1:length(groups)
    group_results{j} = cell(1,3);
    group_results{j}{1} = cell2mat(volumes(groups{j}));
    group_results{j}{2} = cell(1,length(classes));
    group_results{j}{3} = cell(1,length(classes));
    for k = 1:length(classes)
        group_results{j}{2}{1,k} = cell2mat(cellfun( @(a) a(:,k), counts(groups{j}), 'UniformOutput',false));
        group_results{j}{3}{1,k} = cell2mat(cellfun( @(a) a(:,k), densities(groups{j}), 'UniformOutput',false));
    end
end

% statistical analysis
df_stats = get_df_stats(df_template, group_names, classes, group_results);

save(fullfile(result_dir,'stats.mat'), 'df_stats');
if save_table
    save_stats = fullfile(result_dir,'stats.csv');
    writetable(df_stats,save_stats)
end

%% Functions

function [df_volumes, num, density, df_cells] = quantify(data_dir, classes, header, save_table)

    [~, dataset_name] = fileparts(data_dir);

    % Measure Stucture Volumes

    % Take a registered annotation volume and calculate structure voxel volumes. 
    
    resolution = repmat(25,1,3);
    % Calculate volume per voxel in mm^3
    mm = prod(resolution)/(1E3^3);

    % Read structure template file
    df_template = readtable(fullfile(fileparts(which('NM_config')),...
        'annotations','structure_template.csv'));
    indexes = df_template.index;

    % Measure number of voxels for each structure

    mhd_annotation_file = fullfile(data_dir, 'volume\result.mhd');
    [img, ~] = read_mhd(mhd_annotation_file);
    I_mask = permute(img.data,[2,1,3]);
    total_volume = sum(I_mask(:)>0)*mm;

    sums = zeros(1,length(indexes));
    counted = histcounts(I_mask(:),'BinMethod','integers');
    sums(1:length(counted)) = counted;
    df_volumes = sums*mm';

    % Set top row equal to 0 as this is the background
    df_volumes(1) = 0;

    % Sum according to structure level order, except for background
    ids = df_template.id;
    path = df_template.structure_id_path;
    df_new = zeros(size(df_volumes));
    for i = 2:length(ids)
        idx = cellfun(@(s) contains(s,string("/"+ids(i)+"/")), path);
        df_new(i) = sum(df_volumes(idx));
    end
    df_volumes = df_new';

    % Create header name
    df_header = dataset_name + "_Volume";

    % Convert to table
    table_volume = array2table(df_volumes,'VariableNames',df_header);

    % Concatenate counts to results table
    df_vol = horzcat(df_template(:,[1,end-1,end]), table_volume);
    fprintf('%s\t Total structure volume: %.2f mm^3\n',datetime('now'),total_volume)

    % Write new results file - volume

    if save_table
        save_vol = fullfile(data_dir, 'result_volume.csv');
        writetable(df_vol,save_vol)
    end

    % Cell count/density quantification

    num = cell(1,length(classes));
    density = cell(1,length(classes));
    df_cells = df_vol;
    for i = 1:length(classes)
        df_new = zeros(size(df_volumes));
        df_csv = readtable(fullfile(data_dir,'cell_registration',num2str(i),'cells.csv'));
        df_csv.Properties.VariableNames = header;
        count = zeros(1,length(indexes));
        counted = histcounts(df_csv.id,'BinMethod','integers');
        count(1:(length(counted)-1)) = counted(2:end);
        count(1) = 0;
        for j = 2:length(ids)
            idx = cellfun(@(s) contains(s,string("/"+ids(j)+"/")), path);
            count(j) = sum(count(idx));
        end
        df_new(1:length(count)) = count;
        num{1,i} = df_new;
        density{1,i} = num{1,i}./df_volumes;
        table_count = array2table(num{1,i},'VariableNames',classes(i));
        table_density = array2table(density{1,i},'VariableNames',classes(i)+'_density');
        df_cells = horzcat(df_cells,table_count,table_density);
    end
    num = cell2mat(num); density = cell2mat(density);

    % Write new results file - density

    if save_table
        save_cells = fullfile(data_dir,'result_density.csv');
        writetable(df_cells,save_cells)
    end

end

function stats = get_stats(group1,group2,paired,min_stat)

    if nargin<4
        min_stat=0;
    end
    stats = zeros(size(group1,1),8);

    % Mean
    stats(:,1) = mean(group1,2);
    stats(:,2) = mean(group2,2);
    % Standard deviation
    stats(:,3) = std(group1,0,2);
    stats(:,4) = std(group2,0,2);
    % Percent Relative Change
    stats(:,5) = 100*(stats(:,2)-stats(:,1))./stats(:,1);    
    % Set any cells in background to 0
    % stats(1,:) = 0;
    
    % p-value
    if paired
        [~,stats(:,6)] = ttest(group1',group2'); 
    else            
        [~,stats(:,6)] = ttest2(group1',group2','Vartype','unequal'); 
    end
    
    % q-value
    % First remove structures with low counts/no volume
    if nargin<4
        s_pos = ~isnan(stats(:,6));
    else    
        s_pos = max(stats(:,1:2),[],2) > min_stat;
    end
    p_val_thresholded = stats(s_pos,6);
    [stats(s_pos,8), ~, ~, stats(s_pos,7)]=fdr_bh(p_val_thresholded);

    % stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.05);
    stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.01);
    stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.001);
    stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.0001);

    stats(~s_pos,1:5) = 0;
    stats(~s_pos,6:7) = 1;
    stats(~s_pos,8) = 0;

end

function df_stats = get_df_stats(df_template, group_names, classes, group_results)

    df_stats = cell(1,3);
    header_prefix = [repelem(["Mean","StdDev"],2),"PChange","p","p_adj","sig"];
    % Volume
    % Get stats
    stats = get_stats(group_results{1}{1},group_results{2}{1},false,0);
    % Create headers
    df_header = header_prefix + "_" + "Volume";
    df_header(1:4) = repmat(group_names,1,2) + "_" + df_header(1:4);
    % Convert to table
    df_stats{1}  = array2table(stats,'VariableNames',df_header);
    
    % Counts
    % For each channel, calculate stats
    stats = cell(1,length(classes));
    for c = 1:length(classes)
        % Get stats
        sub_stats = get_stats(group_results{1}{2}{c},group_results{2}{2}{c},false,0);
        % Create headers
        df_header = header_prefix + "_" + repelem(classes(c),8) + "_Counts";
        df_header(1:4) = repmat(group_names,1,2) + "_" + df_header(1:4);
        % Convert to table
        stats{c} = array2table(sub_stats,'VariableNames',df_header);
    end
    df_stats{2} = cat(2,stats{:});
    
    % Densities
    % For each channel, calculate stats
    stats = cell(1,length(classes));
    for c = 1:length(classes)
        % Get stats
        sub_stats = get_stats(group_results{1}{3}{c},group_results{2}{3}{c},false,0);
        % Create headers
        df_header = header_prefix + "_" + repelem(classes(c),8) + "_Densities";
        df_header(1:4) = repmat(group_names,1,2) + "_" + df_header(1:4);
        % Convert to table
        stats{c} = array2table(sub_stats,'VariableNames',df_header);
    end
    df_stats{3} = cat(2,stats{:});
    df_stats = horzcat(df_template,cat(2,df_stats{:}));

end
