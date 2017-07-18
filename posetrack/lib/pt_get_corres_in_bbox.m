function [ out_corres, idx ] = pt_get_corres_in_bbox(in_corres, bbox)
%PT_GET_CORRES_AROUND_POINT Summary of this function goes here
%   Detailed explanation goes here

if(isempty(in_corres))
    out_corres = [];
    return;
end

x1 = bbox(1);
y1 = bbox(2);
x2 = bbox(3);
y2 = bbox(4);

idx = bitand(in_corres(:,1) >= x1, in_corres(:,1) <= x2);
tmp_corres = in_corres;
tmp_corres(~idx,:) = nan;

idx = bitand(tmp_corres(:,2) >= y1, tmp_corres(:,2) <= y2);
out_corres = tmp_corres(idx,:);

int_idx = [1:size(in_corres,1)]';
idx = int_idx(idx);

end

