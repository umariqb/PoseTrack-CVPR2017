function pt_track_joint( expidx,cidx,firstidx,nVideos,bRecompute,bVis )
%PT_TRACK_JOINT Summary of this function goes here
%   Detailed explanation goes here

fprintf('pt_track_joint()\n');

if (ischar(expidx))
    expidx = str2num(expidx);
end

if (ischar(cidx))
    cidx = str2num(cidx);
end

if (ischar(firstidx))
    firstidx = str2num(firstidx);
end

if (nargin < 2)
    firstidx = 1;
end

if (nargin < 3)
    nVideos = 1;
elseif ischar(nVideos)
    nVideos = str2num(nVideos);
end

if (nargin < 4)
    bRecompute = false;
end

if (nargin < 5)
    bVis = false;
end

fprintf('expidx: %d\n',expidx);
fprintf('cidx: %d\n', cidx);
fprintf('firstidx: %d\n',firstidx);
fprintf('nImgs: %d\n',nVideos);

is_release_mode = get_release_mode();

video_set = 'test';
p = pt_exp_params(expidx);
exp_dir = fullfile(p.expDir, p.shortName);
load(p.testGT,'annolist');

ptMulticutDir = p.ptMulticutDir;
mkdir_if_missing(ptMulticutDir);
fprintf('multicutDir: %s\n',ptMulticutDir);
visDir = fullfile(ptMulticutDir, 'vis');
mkdir_if_missing(visDir);

num_videos = length(annolist);

lastidx = firstidx + nVideos - 1;
if (lastidx > num_videos)
    lastidx = num_videos;
end

if (firstidx > lastidx)
    return;
end

% computation parameters
pairwiseDir = p.pairwiseDir;
ptPairwiseDir = p.ptPairwiseDir;
pad_orig = p.([video_set 'Pad']);


fprintf('Loading temporal model from %s\n', ptPairwiseDir);
modelName  = [ptPairwiseDir '/temporal_model_cidx_' num2str(cidx)];
temporal_model = struct;
m = load(modelName,'temporal_model');
temporal_model.diff = struct;
temporal_model.diff.log_reg = m.temporal_model.log_reg;
temporal_model.diff.training_opts = m.temporal_model.training_opts;

fprintf('recompute %d\n', bRecompute);

idxStart = 1;
for vIdx = firstidx:lastidx
    fprintf('vididx: %d\n',vIdx);
    ann = annolist(vIdx);
    vid_name    = ann.name;
    vid_dir     = fullfile(p.vidDir, vid_name);
    num_frames  = ann.num_frames;
    fn          = dir([vid_dir,'/*.jpg']);
    frame_pairs = pt_build_frame_pairs(num_frames, p.maxFrameDist);
    corres_dir = fullfile(p.correspondences, vid_name);
    scoremaps_dir = fullfile(p.scoremaps, vid_name);
    
    detPerFrame = zeros(num_frames,1); 
    detections = [];
    for fIdx = 1:num_frames
        fprintf('Frame: %d/%d\n',fIdx, num_frames);
        fr_fn = fullfile(vid_dir, fn(fIdx).name);
        [~,fr_name,~] = fileparts(fr_fn);        
        score_fn = fullfile(scoremaps_dir, fr_name);

        orig_frame = imread(fr_fn);
        [w,h,~] = size(orig_frame);        
        max_side = max(h,w);
        
        if(max_side > p.maxDim)
            scale = p.maxDim/max_side;
            frame = imresize(orig_frame, scale);
        else
            frame = orig_frame;
        end
        
        load(score_fn, 'scoremaps');
        
        if bVis && false
            figure(12);
            for s=1:length(p.scales)
                figure(12);
                smaps = scoremaps{s};
                sm = smaps(:,:,cidx);
                imagesc(sm);
                colorbar;
                pause();
            end
        end
        
        dets = pt_scoremap_to_detections(p, score_fn, cidx);
        num_dets = size(dets.unProb, 1);
        dets.frameIndex = ones(num_dets,1)*fIdx;
        dets.index = [idxStart:idxStart+num_dets-1]';
        idxStart = idxStart + num_dets;
        
        if (bVis)
            figure(100);clf;
            imagesc(frame); axis equal; hold on;
            %[val,idx] = sort(dets.unProb(:,j),'descend');
            cb = dets.unPos;
            scale = dets.scale;
            plot(cb(:,1),cb(:,2),'b+');
            for j = 1:size(dets.unPos, 1)
                text(cb(j,1),cb(j,2), num2str(j), 'Color', 'g', 'FontSize', 10);
                r = (p.patchSize * 1/scale(j))/2;
                rectangle('Position', [cb(j,1)-r,cb(j,2)-r, 2*r, 2*r]);
            end
            pause(0.01);

            if ~is_release_mode
                pause;
            end
        end
        detPerFrame(fIdx) = size(dets.unProb,1);
        
        if(isempty(detections))
            detections = dets;
        else            
            detections = Detections.merge(detections, dets);
        end
    end
      
    % compute pairwise probabilties
    pwProb = pt_compute_temporal_pairwise_probabilities(p, detections, temporal_model, cidx, fn, corres_dir);
    
    pwProbSolver = pwProb(:,3:6);
    pwProbSolver(:,1:2) = pwProbSolver(:,1:2) - 1;
    numDetections = size(detections.unProb, 1);
        
    problemFname = [ptMulticutDir '/problem-' vid_name '.h5'];
    solutionFname = [ptMulticutDir '/solution-' vid_name '.h5'];

    fprintf('save problem\n');
    % write problem
    dataName    = 'detections-info';
    write_mode  = 'overwrite';
    marray_save(problemFname, dataName, numDetections, write_mode);
    
    dataName = 'pt-join-probabilities';
    write_mode = 'append';
    marray_save(problemFname, dataName, pwProbSolver, write_mode);

    solver = p.ptSolver;
    time_limit = p.time_limit;
    cmd = [solver ' ' problemFname  '  ' solutionFname ' ' num2str(time_limit)];
    
    fprintf('calling pt-solver: %s\n', cmd);

    [~,hostname] = unix('echo $HOSTNAME');
            
    fprintf('hostname: %s',hostname);
    pre_cmd = ['export GRB_LICENSE_FILE=' p.gurobi_license_file];

    tic
    setenv('LD_LIBRARY_PATH', '');
    s = unix([pre_cmd '; ' cmd]);
    toc
    if (s > 0)
        error('solver error');
    end
    assert(s == 0);

    % clean up
    unix(['rm ' problemFname]);

    % load solution
    dataName = 'part-tracks';
    unLab = marray_load(solutionFname, dataName);
    
    pt_visualize_tracks(p, detections, unLab+1, vid_dir, num_frames);

end

end

