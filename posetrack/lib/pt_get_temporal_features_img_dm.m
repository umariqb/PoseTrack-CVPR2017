function [ feat ] = pt_get_temporal_features_img_dm(p, locs1, corres_pts1, locs2, corres_pts2, idxs_pt_pair )
%PT_GET_TEMPORAL_FEATURES_IMG_DM Summary of this function goes here
%   Detailed explanation goes here


if (nargin < 5)
    [p,q] = meshgrid(size(locs1,1), size(locs2,1));
    idxs_pt_pair = [p(:) q(:)];
end

if(isempty(locs1) || isempty(locs2) || isempty(corres_pts1) || isempty(corres_pts2))
    feat = [];
    return;
end

num_pairs = size(idxs_pt_pair,1);

feat =  zeros(num_pairs, 11);

try
    
for i = 1:num_pairs
    
    pt1       = locs1(idxs_pt_pair(i,1), 1:2);
    scale1    = 1/locs1(idxs_pt_pair(i,1), 3);
    score1    = locs1(idxs_pt_pair(i,1), end);
    [~, idx1] = pt_get_corres_around_point(corres_pts1, pt1, p.patchSize*scale1);
    
    pt2       = locs2(idxs_pt_pair(i,2), 1:2);
    scale2    = 1/locs2(idxs_pt_pair(i,2), 3);
    score2    = locs2(idxs_pt_pair(i,2), end);
    [~, idx2] = pt_get_corres_around_point(corres_pts2, pt2, p.patchSize*scale2);
    
    inter = intersect(idx1, idx2);
    uni   = union(idx1, idx2);
    
    if(~isempty(idx1) || ~isempty(idx2))
        feat(i, 1) = length(inter)/length(uni);
    end
    feat(i, 2) = min(score1, score2);
    feat(i, 3) = feat(i,1)*feat(i,2);
    feat(i, 4) = feat(i,1).^2;
    feat(i, 5) = feat(i,2).^2;
    
    dist = pt1*scale1 - pt2*scale1; % multiply only with one of the scales.
    feat(i,6:7) = dist;
    feat(i,8) = norm(dist);
    feat(i,9:10) = feat(i,6:7).^2;
    feat(i,11) = feat(i,8)^2;
    
end

catch
    keyboard();
end
end

