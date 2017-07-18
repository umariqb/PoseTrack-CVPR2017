function pt_demo_multiperson()

experiment_index = 1;
video_index = 1;
video_set = 'test';

set_release_mode(true);


% process input image with the CNN and cache confidence maps to the disk
pt_cache_cnn_features(experiment_index, video_set, video_index, 22);

% extract deep-matching features and cache to the disk.
pt_cache_deepmatching_features(experiment_index, video_set, video_index, 22);

% prepare and run ILP inference
test_spatial_app_neighbour(experiment_index, image_index, 1, true, true);

% visualise predictions
vis_people(experiment_index, image_index);

end



