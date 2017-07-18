function pt_pose_tracking_nodes_and_edges( expidx,firstidx,nVideos,bRecompute,bVis)
%PT_TRACK_JOINT Summary of this function goes here
%   Detailed explanation goes here

fprintf('pt_pose_tracking()\n');

if (ischar(expidx))
    expidx = str2num(expidx);
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

ptDetectionDir = p.ptDetectionsDir;
mkdir_if_missing(ptDetectionDir);
fprintf('detectionDir: %s\n',ptMulticutDir);

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
stride = p.stride;
half_stride = stride/2;
locref_scale_mul = p.locref_scale;
nextreg = isfield(p, 'nextreg') && p.nextreg;
unLab_cls = 'uint64';
max_sizet_val = intmax(unLab_cls);
ptPairwiseDir = p.ptPairwiseDir;
pad_orig = p.([video_set 'Pad']);

temporal_models = cell(14,1); % head for now
temporal_cidxs = p.temporal_cidxs;
for cidx=temporal_cidxs
    fprintf('Loading temporal model from %s\n', ptPairwiseDir);
    modelName  = [ptPairwiseDir '/temporal_model_with_dist_cidx_' num2str(cidx)];
    temporal_model = struct;
    m = load(modelName,'temporal_model');
    temporal_model.diff = struct;
    temporal_model.diff.log_reg = m.temporal_model.log_reg;
    temporal_model.diff.training_opts = m.temporal_model.training_opts;
    temporal_models{cidx}=temporal_model;
end

stagewise = false;
if isfield(p, 'cidxs_full')
    cidxs_full = p.cidxs_full;
else
    cidxs_full = p.cidxs;
end

if stagewise
    num_stages = length(p.cidxs_stages);
else
    num_stages = 1;
end

if isfield(p, 'dets_per_part_per_stage')
    dets_per_part = p.dets_per_part_per_stage;
elseif isfield(p, 'dets_per_part')
    dets_per_part = p.dets_per_part;
    if stagewise
        dets_per_part = repmat(dets_per_part, num_stages);
    end
end

use_graph_subset = isfield(p, 'graph_subset');
if use_graph_subset
    load(p.graph_subset);
end

pairwise = load_pairwise_data(p);
graph = pairwise.graph;

if isfield(p, 'nms_dist')
    nms_dist = p.nms_dist;
else
    nms_dist = 1.5;
end

hist_pairwise = isfield(p, 'histogram_pairwise') && p.histogram_pairwise;

pwIdxsAllrel1 = build_pairwise_pairs(cidxs_full);

fprintf('Loading spatial model from %s\n', pairwiseDir);
for sidx1 = 1:length(pwIdxsAllrel1)
    cidx1 = pwIdxsAllrel1{sidx1}(1);
    cidx2 = pwIdxsAllrel1{sidx1}(2);
    modelName  = [pairwiseDir '/spatial_model_cidx_' num2str(cidx1) '_' num2str(cidx2)];
    spatial_model{sidx1} = struct;
    if hist_pairwise
        spatial_model{sidx1}.diff = load(modelName,'edges1','edges2','posHist','negHist','nPos','nNeg');
    else
        m = load(modelName,'spatial_model');
        spatial_model{sidx1}.diff = struct;
        spatial_model{sidx1}.diff.log_reg = m.spatial_model.log_reg;
        spatial_model{sidx1}.diff.training_opts = m.spatial_model.training_opts;
    end
end

fprintf('recompute %d\n', bRecompute);

nNodes = 0;
nSpatialEdges = 0;
nTemporalEdges = 0;

bVis = true;
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
    
    cidxs = cidxs_full;
        
    fname = [ptMulticutDir '/' vid_name '_cidx_' num2str(cidxs(1)) '_' num2str(cidxs(end))];
    detFname = [p.ptDetectionsDir '/' vid_name];
    
    predFname = [ptMulticutDir '/prediction_' num2str(expidx) '_' vid_name '.mat'];
    [pathstr,~,~] = fileparts(predFname);
    mkdir_if_missing(pathstr);

%     if(exist(predFname, 'file'))
%         fprintf('Prediction already exist at: %s\n', predFname);
%         continue;
%     end
        
    if(exist([detFname,'.mat'], 'file') && ~bRecompute)
        load(detFname, 'detections');
    else
        [pathstr,~,~] = fileparts(detFname);
        mkdir_if_missing(pathstr);
        for fIdx = 1:50 %num_frames
            fprintf('Frame: %d/%d\n',fIdx, num_frames);
            fr_fn = fullfile(vid_dir, fn(fIdx).name);
            [~,fr_name,~] = fileparts(fr_fn);        

            score_fn = fullfile(scoremaps_dir, fr_name);
            load(score_fn, 'scoremaps');

            frame = imread(fr_fn);
            [w,h,~] = size(frame);        
            max_side = max(h,w);

%             if(max_side > p.maxDim)
%                 scale = p.maxDim/max_side;
%             end

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

            dets = pt_scoremap_to_detections(p, score_fn, cidxs);
            num_dets = size(dets.unProb, 1);
            dets.frameIndex = ones(num_dets,1)*fIdx;
            dets.index = [idxStart:idxStart+num_dets-1]';
            idxStart = idxStart + num_dets;
        
            if (bVis && false)
                figure(100);clf;
                imagesc(frame); axis equal; hold on;
                %[val,idx] = sort(dets.unProb(:,j),'descend');
                cb = dets.unPos;
                scale = dets.scale;
                plot(cb(:,1),cb(:,2),'b+');
%                 for j = 1:size(dets.unPos, 1)
% %                     text(cb(j,1),cb(j,2), num2str(j), 'Color', 'g', 'FontSize', 10);
%                     r = (p.patchSize * 1/scale(j))/2;
% %                         rectangle('Position', [cb(j,1)-r,cb(j,2)-r, 2*r, 2*r]);
%                 end
                pause(0.01);

%                 if ~is_release_mode
% %                     pause;
%                 end
            end
            detPerFrame(fIdx) = size(dets.unProb,1);

            if(isempty(detections))
                detections = dets;
            else            
                detections = Detections.merge(detections, dets);
            end
        end
            save(detFname, 'detections');
    end
        
    idxs = ismember(detections.partClass, cidxs);
    detections = Detections.slice(detections, idxs);
    
    temporalWindows = pt_generate_frame_windows(num_frames, p.temporalWinSize);
    
    for w = 1:1%size(temporalWindows, 1)
        stIdx  = temporalWindows(w, 1);
        endIdx = temporalWindows(w, 2);
            
        idxs = ismember(detections.frameIndex, stIdx:endIdx) > 0;
%         idxs1 = bitand(detections.unPos(:,1) > 303, detections.unPos(:,1) < 542);
%         idxs2 = bitand(detections.unPos(:,2) > 148, detections.unPos(:,2) < 430);
%         idxs = bitand(idxs, idxs1);
%         idxs = bitand(idxs, idxs2);
        wDetections = Detections.slice(detections, idxs);
        wDetections.unLab = zeros(size(wDetections.unProb,1),2,unLab_cls);
        wDetections.unLab(:,:) = max_sizet_val;

        if(w > 1)
            fprintf('Window %d/%d: Previous: %d \t New: %d\n', w, size(temporalWindows, 1), size(prev_dets.unProb,1), size(wDetections.unProb,1));
            wDetections = Detections.merge(prev_dets, wDetections); 
        end
        wDetections.index = [1:size(wDetections.unProb,1)]';
        fprintf('Number of Detections: %d\n', size(wDetections.unProb,1));
                    
        pwProbTemporal = [];
        for cidx=temporal_cidxs
            idx = wDetections.partClass == cidx;
%             if(w > 1)
%                 % compute temporal features only for the overlap. 
%                 idx = wDetections.partClass == 14;
%                 oIdx = ismember(wDetections.frameIndex, max(1,stIdx-p.maxFrameDist:endIdx)) > 0;
%                 idx = bitand(idx,oIdx);
%             end
            fprintf('Computing temporal pairwise probabilities. Part = %d\n',cidx);
            tempDetections = Detections.slice(wDetections, idx);
            
             pwProb = pt_compute_temporal_pairwise_probabilities(p, tempDetections, temporal_models{cidx}, cidx, fn, corres_dir);
%             if(ismember(cidx, [9,10,3,4]))
%                 pwProb(:,end) = pwProb(:,end).^3;
%             end

            pwProbTemporal = [pwProbTemporal;pwProb];
        end
        % compute spatial pairwise probabilities
        fprintf('Computing spatial pairwise probabilities.\n');
        pwProbSpatial = cell(num_frames,1);
        for fIdx = 1:endIdx
            idx = wDetections.frameIndex == fIdx;
            fDetections = Detections.slice(wDetections, idx);
            if(size(fDetections.unProb,1) == 0)
                continue;
            end
            if (bVis && false)
                fr_fn = fullfile(vid_dir, fn(fIdx).name);
                frame = imread(fr_fn);
                figure(101+fIdx);clf;
                imagesc(frame); axis equal; hold on;
                %[val,idx] = sort(dets.unProb(:,j),'descend');
                cb = fDetections.unPos;
                scale = fDetections.scale;
%                 plot(cb(:,1),cb(:,2),'b+');
                for j = 1:size(fDetections.unPos, 1)
                    if(ismember(fDetections.partClass(j),[9,13]))
                        text(cb(j,1),cb(j,2), num2str(fDetections.index(j)), 'Color', 'g', 'FontSize', 10);
                        r = (p.patchSize * 1/scale(j))/2;
                        rectangle('Position', [cb(j,1)-r,cb(j,2)-r, 2*r, 2*r]);
                    end
                end
                pause();
            end

            pwProbSpatial{fIdx,1} = ... 
                pt_compute_spatial_pairwise_probabilities( p, fDetections,  ... 
                    cidxs, pwIdxsAllrel1, graph, hist_pairwise, nextreg, spatial_model, ...
                    use_graph_subset);

            pwProbSpatial{fIdx,1} = cat(2,ones(size( pwProbSpatial{fIdx,1}, 1),1)*fIdx, pwProbSpatial{fIdx,1});
        end
        pwProbSpatial = cell2mat(pwProbSpatial);

        % prepare problem for solver
        numDetections = size(wDetections.unProb, 1);
        pwProbSolverTemp = pwProbTemporal(:,3:6);
        pwProbSolverTemp(:,1:2) = pwProbSolverTemp(:,1:2) - min(wDetections.index);
%         idxs = pwProbSolverTemp(:,3) <  0.6;
        pwProbSolverTemp(:,3) = pwProbSolverTemp(:,3);
        pwProbSolverSpat = pwProbSpatial(:, 2:6);
        pwProbSolverSpat(:,1:2) = pwProbSolverSpat(:,1:2) - min(wDetections.index); 

        
        nNodes(vIdx)          =  numDetections;
        nSpatialEdges(vIdx)   =  size(pwProbSolverSpat,1);
        nTemporalEdges(vIdx)  =  size(pwProbSolverTemp,1); 
    end 
end  

fprintf('Average Number of Nodes = %10.0d\n', uint16(mean(nNodes)));
fprintf('Average Number of Spatial Edges = %d\n', uint16(mean(nSpatialEdges)));
fprintf('Average Number of Temporal Edges = %d\n', uint16(mean(nTemporalEdges)));
fprintf('Median Number of Nodes = %d\n', uint16(median(nNodes)));
fprintf('Median Number of Spatial Edges = %d\n', uint16(median(nSpatialEdges)));
fprintf('Median Number of Temporal = %d\n', uint16(median(nTemporalEdges)));




fprintf('done\n');

if (isdeployed)
    close all;
end
end


function [outClusters] = pruneClusters(clusters, frameIndexes, labels, minLen)
    outClusters = [];
    for c = 1:length(clusters)
        cl = clusters(c);
        idxs = labels == cl;
        frameLen = length(unique(frameIndexes(idxs)));
        if(frameLen > minLen)
            outClusters = [outClusters;cl];
        end
    end
end

function [unLab_new, idxs] = pruneDetections(unLab, detections, cidxs, minLen)
    
    clusters = unique(unLab(:, 2));
    clusters = pruneClusters(clusters, detections.frameIndex, unLab(:,2), minLen);
    idxs = [];
    
    stIdx = min(detections.frameIndex);
    endIdx = max(detections.frameIndex);
    
    for fIdx = stIdx:endIdx
        
        %get the detections in current frame and suppress those with status 0
        cfIdxs = find(bitand(detections.frameIndex == fIdx, logical(unLab(:,1))));         
        
        for j = 1:length(clusters)
            cl = clusters(j);
            cDets  = unLab(cfIdxs, 2) == cl;
            labels = detections.partClass(cfIdxs);
            for k = 1:length(cidxs)
                cidx = cidxs(k);
                bundle = find(cDets & labels == cidx);
                if isempty(bundle)
                    continue;
                end
                probs = detections.unProb(cfIdxs(bundle));
                [~, I] = max(probs);
                idxs = [idxs; cfIdxs(bundle(I))];
            end
        end 
    end
    
    idxs = sort(idxs);
    
    unLab_new = unLab(idxs, :);
    for j = 1:length(clusters)
        cl = clusters(j);
        I = unLab_new == cl;
        unLab_new(I) = j-1;
    end
        
    %idxs
end

function dets = copy_detections(dets_src)
    dets = struct();
    dets.unPos = dets_src.unPos;
    dets.unPos_sm = dets_src.unPos_sm;
    dets.unProb = dets_src.unProb;
    dets.locationRefine = dets_src.locationRefine;
    dets.nextReg = dets_src.nextReg;
    dets.unProbNoThresh = dets_src.unProbNoThresh;
    if isfield(dets_src, 'unLab')
        dets.unLab = dets_src.unLab;
    end
    if isfield(dets_src, 'scale')
        dets.scale = dets_src.scale;
    end
    if isfield(dets_src, 'frameIndex')
        dets.frameIndex = dets_src.frameIndex;
    end
    if isfield(dets_src, 'index')
        dets.index = dets_src.index;
    end
    if isfield(dets_src, 'partClass')
        dets.partClass = dets_src.partClass;
    end

end


