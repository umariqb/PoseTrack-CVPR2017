function pt_cache_cnn_features( expidx, video_set, firstidx, nVids, bRecompute)

p = pt_exp_params(expidx);

if (nargin < 3)
    firstidx = 1;
elseif ischar(firstidx)
    firstidx = str2num(firstidx);
end

if(nargin < 5)
    bRecompute = false;
end

if strcmp(video_set, 'test')
    load(p.testGT)
else
    load(p.trainGT)
end

num_videos = length(annolist);

if (nargin < 4)
    nVids = num_videos;
elseif ischar(nVids)
    nVids = str2num(nVids);
end

lastidx = firstidx + nVids - 1;
if (lastidx > num_videos)
    lastidx = num_videos;
end


mirror_map = [6 5 4 3 2 1 12 11 10 9 8 7 13 14];

nextreg = isfield(p, 'nextreg') && p.nextreg;

models_dir = [p.code_dir '/deepcut/models/'];
net_def_file = [models_dir p.net_def_file];

net_dir = p.net_dir;
net_bin_file = get_net_filename(net_dir);
if isfield(p, 'net_bin_file')
    net_bin_file = [net_dir '/' p.net_bin_file];
end

caffe.reset_all();
caffe.set_mode_gpu();
net = caffe.Net(net_def_file, net_bin_file, 'test');
fprintf('testing from net file %s\n', net_bin_file);

pairwise = load_pairwise_data(p);

for vid_idx = firstidx:lastidx

    vid_name = annolist(vid_idx).name;
    
    cache_dir = fullfile(p.scoremaps, vid_name);
    mkdir_if_missing(cache_dir);
    fprintf('save dir %s\n', cache_dir);
    
    vid_dir = fullfile(p.vidDir, vid_name);
    num_frames = annolist(vid_idx).num_frames;
    fn = dir([vid_dir,'/*.jpg']);
    
    if(length(fn) ~= num_frames)
        keyboard();
    end
    
    for fr_idx =1:num_frames
        fprintf('%s: test (%s) %d/%d\n', procid(), p.name, fr_idx, num_frames);
        fr_fn = fullfile(vid_dir, fn(fr_idx).name);
        [~,fr_name,~] = fileparts(fr_fn);
        save_file_name = fullfile(cache_dir, fr_name);
    
        if (exist([save_file_name '.mat'], 'file') == 2) && ~bRecompute
            continue
        end

        im = imread(fr_fn);

        pad_orig = p.([video_set 'Pad']);
        [w,h,~] = size(im);
        
        max_side = max(h,w);
        
        if(max_side > p.maxDim)
            scale = p.maxDim/max_side;
            im = imresize(im, scale);
            [w,h,~] = size(im);
        end
        
        crop = [1,1,h,w];        
        if(p.multiscale)
            [scoremaps, locreg_pred, nextreg_pred] = pt_extract_features_multiscale(im, net, p, [], pad_orig, pairwise, crop);
            if(false);
                for s=1:length(p.scales)
                    figure(13);
                    imshow(im);
                    figure(12);
                    scmap_vis = visualise_scoremap(scoremaps{s});
                    imshow(scmap_vis);
                    pause();
                end
            end
        else
            [scoremaps, locreg_pred, nextreg_pred] = extract_features(im, net, p, [], pad_orig, pairwise, crop);
            if(false);
                figure(12);
                scmap_vis = visualise_scoremap(scoremaps);
                imshow(scmap_vis);
                pause();
            end
        end
        save(save_file_name, 'scoremaps');
        if ~isempty(locreg_pred)
            save(fullfile(cache_dir, [fr_name '_locreg']), 'locreg_pred');
        end
        if nextreg && ~isempty(nextreg_pred)
            save(fullfile(cache_dir, [fr_name '_nextreg']), 'nextreg_pred');
        end
    end
end

caffe.reset_all();
end