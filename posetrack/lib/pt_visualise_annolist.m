function pt_visualise_annolist( expidx, image_set, firstidx )

if (nargin < 3)
    firstidx = 1;
end

% index: frame number
% person_ids: ids of the persons to be plotted
cmap = hsv(14);  
lineWidth = 3;
edges = [1 2; 2 3; 4 5; 5 6; 4 7; 3 7; 7 8; 9 10; 10 11; 7 11; 7 12; 12 13; 13 14];

p = pt_exp_params(expidx);

load(p.([image_set 'GT']), 'annolist');

num_videos = length(annolist);
parts = get_parts();

for vididx = firstidx:num_videos
    fprintf('Video %d/%d\n', vididx, num_videos);
    
    vid_dir = fullfile(p.vidDir, annolist(vididx).name);
    num_frames = annolist(vididx).num_frames;
    fn = dir([vid_dir,'/*.jpg']);
    
    for imgidx = 1:num_frames
        
        img_fn = fullfile(vid_dir, fn(imgidx).name);
    
        im = imread(img_fn);
        clf;
        imagesc(im); axis equal; hold on;
    
        for k = 1:annolist(vididx).num_persons
            joints = pt_get_anno_joints(annolist(vididx).annopoints{k,imgidx}, p.pidxs, parts);
            figure(1);
            vis_pred(joints);
            if(isfield(annolist(vididx).annopoints{k,imgidx}, 'head_rect'))
                rect = annolist(vididx).annopoints{k,imgidx}.head_rect;
                rectangle('Position', [rect(1),rect(2),rect(3)-rect(1),rect(4)-rect(2)], 'EdgeColor', 'r');
            end
        end
        pause(0.001);
        hold off;
    %{
%     [scmap, poly] = get_sticks_segmentation(p, im, joints);
%     plot(poly(:,1), poly(:,2));
%     
%     scmap = visualise_scoremap(scmap , 4);
%     figure(2);
%     imshow(scmap);
    %}

    end
end


end

