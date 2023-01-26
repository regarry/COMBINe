% perform channel alignment for a single dataset

sample = 'sample name as in NM_samples';
[img_directory, output_directory, group] = NM_samples(sample, 'True');
[config, path_table] = NM_config('process', sample); % config dataset
[config, path_table] = perform_channel_alignment(config, path_table, 'True'); % channel alignment
translation_table = reformat(config); % reformat output for Python loading