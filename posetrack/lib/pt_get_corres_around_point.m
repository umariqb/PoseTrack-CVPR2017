function [ out_corres, idx ] = pt_get_corres_around_point(in_corres, pt, patch_size)
%PT_GET_CORRES_AROUND_POINT Summary of this function goes here
%   Detailed explanation goes here

if(isempty(in_corres))
    out_corres = [];
    return;
end

r = patch_size / 2;

idx = bitand(in_corres(:,1) >= pt(1)-r, in_corres(:,1) <= pt(1)+r);
tmp_corres = in_corres;
tmp_corres(~idx,:) = nan;

idx = bitand(tmp_corres(:,2) >= pt(2)-r, tmp_corres(:,2) <= pt(2)+r);
out_corres = tmp_corres(idx,:);

int_idx = [1:size(in_corres,1)]';
idx = int_idx(idx);

end

