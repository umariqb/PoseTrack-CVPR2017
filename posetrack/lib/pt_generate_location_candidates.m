function [candidates ] = pt_generate_location_candidates(filename, p, cidx)

load(filename, 'scoremaps')

assert(p.multiscale && iscell(scoremaps) && ... 
        length(scoremaps) == length(p.scales));


num_scales = length(p.scales);
stride = p.stride;

locations = [];
scores    = [];
scales    = [];

for s = 1:num_scales
    
    cur_scale = p.scales(s);
    inv_cur_scale = 1/cur_scale;
    
    s_scoremaps = scoremaps{s};
    height = size(s_scoremaps, 1);
    width  = size(s_scoremaps, 2);
    
    num_candidates = height*width;
    num_scores = size(s_scoremaps, 3);
    size_2d = size(s_scoremaps(:,:,1));  
    
    s_locations = zeros(num_candidates, 2);
    s_scores = zeros(num_candidates, length(cidx));
    s_scales = zeros(num_candidates, 1);
    for jj = 1:height
        for ii = 1:width
            ind = sub2ind(size_2d, jj, ii);
            s_locations(ind, :) = [ii-1, jj-1]*stride*inv_cur_scale;
            s_scores(ind, :) = s_scoremaps(jj, ii, cidx);
            s_scales(ind, :) = cur_scale;
        end
    end
    
    locations = [locations; s_locations];
    scores    = [scores; s_scores];
    scales  = [scales; s_scales];    
end

candidates = cat(2, locations, scales, scores);
end


