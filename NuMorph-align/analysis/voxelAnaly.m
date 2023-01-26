% Voxel-wise quantitative analysis between groups

%% Create voxel image for each class

dataset_dir = 'path\to\all\datasets\and\quantifications';
classes=["yellow_neuron","yellow_astrocyte","green_neuron",...
    "green_astrocyte","red_neuron","red_astrocyte"];

result_dir = fullfile(dataset_dir, 'analysis');
mkdir(result_dir);
[data_subdirs,~] = GetSubDirsFirstLevelOnly(dataset_dir);

% voxel-wise stats

fdr_thresh = 0.05;
orient = [[2,1,3],[2,1,3],[2,1,3],[2,1,3],[-2,1,3],[-2,1,3]]; % same as orientation used in ClearMap registration
voxel_volume = cell(1,length(data_subdirs));
options.overwrite = true;

for i = 1:length(data_subdirs)
    data_dir = fullfile(data_subdirs{i},'cell_registration');
    img = cell(1,length(classes));
    for j = 1:length(classes)
        img_counts = loadtiff(fullfile(data_dir,num2str(j),'density_counts.tif'));
        img_counts = ipermute_orient(img_counts,orient(i));
        img_size = size(img_counts)/4;
        % downsample by 4 (100 um)
        sumrows = sum(reshape(img_counts,4,[]),1); 
        sumcols = sum(reshape(sumrows,img_size(1),4,[]),2);    
        img_counts = reshape(sum(reshape(sumcols,img_size(1)*img_size(2),4,[]),2),img_size);
        saveastiff(single(img_counts),fullfile(data_dir,num2str(j),'density_counts_downsampled.tif'), options);
        img{j} = img_counts;
    end
    voxel_volume{i} = cellfun(@(s) s(:),img,'UniformOutput',false);
end

voxel_volume = cellfun(@(s) cell2mat(s), voxel_volume, 'UniformOutput', false);

groups = {[1,4,6],[2,3,5]}; % assign group indexes and divide into two groups
vox_1 = voxel_volume(groups{1});
vox_2 = voxel_volume(groups{2});

save(fullfile(result_dir,'voxel_vol_100um.mat'), 'voxel_volume');

for i = 1:length(classes)
    
    vox_res = ones(prod(img_size),3);

    % Seperate groups
    set1 = cat(2,cell2mat(cellfun(@(s) s(:,i),vox_1,'UniformOutput',false)));
    set2 = cat(2,cell2mat(cellfun(@(s) s(:,i),vox_2,'UniformOutput',false)));

    % Calculate means
    m_set1 = mean(set1,2);
    m_set2 = mean(set2,2);
    
    idx = find(m_set1>0 & m_set2>0);
    
    % Fold change
    vox_res(:,1) = 0;
    vox_res(idx,1) = (m_set2(idx)-m_set1(idx))./m_set1(idx); 

    % p value
    [~,p] = ttest2(set1(idx,:)',set2(idx,:)');
    vox_res(idx,2) = -log10(p);

    % Adjusted p value. Apply threshold to adjusted p value
    [~,~,~, adj_p] = fdr_bh(p);
    adj_p = -log10(adj_p);
    adj_p(adj_p<-log10(fdr_thresh)) = -log10(fdr_thresh);
    vox_res(idx,3) = adj_p;
    
    % Reshape and save to file
    vox_res = reshape(vox_res,[img_size,3]);
    
    % Save file
    fname = fullfile(result_dir,sprintf("%d_voxels.nii", i));
    niftiwrite(vox_res,fname)
    
end
