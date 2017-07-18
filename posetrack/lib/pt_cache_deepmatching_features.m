function pt_cache_deepmatching_features(expidx, video_set, firstidx, nVids)

p = pt_exp_params(expidx);

if (nargin < 3)
    firstidx = 1;
elseif ischar(firstidx)
    firstidx = str2num(firstidx);
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

% params
overwrite = false;

for vid_idx = firstidx:lastidx

    vid_name = annolist(vid_idx).name;
    
    cache_dir = fullfile(p.correspondences, vid_name);
    mkdir_if_missing(cache_dir);
    fprintf('save dir %s\n', cache_dir);
    
    vid_dir = fullfile(p.vidDir, vid_name);
    num_frames = annolist(vid_idx).num_frames;
    fn = dir([vid_dir,'/*.jpg']);
    
    assert(length(fn) == num_frames)

    frame_pairs = pt_build_frame_pairs(num_frames, p.maxFrameDist);
    
    deepmatching = p.deepMatching;
    
    num_pairs = length(frame_pairs);
    for idx =1:num_pairs
        
        
        pair = frame_pairs(idx,:);
        
        fr_fn1 = fullfile(vid_dir, fn(pair(1)).name);
        fr_fn2 = fullfile(vid_dir, fn(pair(2)).name);
        
        [~,fr_name1,~] = fileparts(fr_fn1);
        [~,fr_name2,~] = fileparts(fr_fn2);
        
        save_file_name = fullfile(cache_dir, [fr_name1,'_',fr_name2,'.txt']);
            
        if (exist([save_file_name], 'file') == 2) && ~overwrite
            continue
        end

        time_start = tic;
        dist = abs(pair(1) - pair(2));
%         rad = p.matchRadius * dist;
        im_sz = size(imread(fr_fn1));
        max_side = max(im_sz);
        
        sf = 1;
        if(max_side > p.maxDimDM)
            sf = 2;
        end
                
        cmd = sprintf('%s %s %s -nt 0 -downscale %d > %s', deepmatching, fr_fn1, fr_fn2, sf, save_file_name);
              
        s = unix(cmd);
        time_elasped = toc(time_start);
        if (s > 0)
            error('DeepMatching error');
        end
        assert(s == 0);
        fprintf('%s: DeepMatching (%s) %d/%d \t elasped time: %f\n', procid(), p.name, idx, num_pairs, time_elasped);
    end
end

end