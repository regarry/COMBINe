%% batch processing

sample = {'sample 1', 'sample 2', 'sample n'}; % names as in NM_samples

for i = 1:length(samples)
    sample = samples{i};
    [img_directory, output_directory, group] = NM_samples(sample, 'True');
    [config, path_table] = NM_config('process', sample); % config dataset
    [config, path_table] = perform_channel_alignment(config, path_table, 'True'); % channel alignment
    translation_table = reformat(config); % reformat output for Python loading
end 