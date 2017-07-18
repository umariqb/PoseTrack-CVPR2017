function [ idxsRel ] = pt_build_frame_pairs( num_frames, max_dist )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[q1,q2] = meshgrid(1:num_frames,1:num_frames);
idxsAllrel = [q1(:) q2(:)];

idxsExc = idxsAllrel(:,1) >= idxsAllrel(:,2);
idxsRel = idxsAllrel(~idxsExc,:);

idxsExc = abs(idxsRel(:,1) - idxsRel(:,2)) > max_dist;
idxsRel = idxsRel(~idxsExc,:);


% cidxs = 1:num_frames;
% pwIdxsAllrel1 = cell(0);
% n = 0;
% for i = 1:length(cidxs)-1
%   for j = i+1:length(cidxs)
%     dist = abs(cidxs(i) - cidxs(j));
%     if(dist > max_dist)
%         continue;
%     end
%     n = n + 1;
%     pwIdxsAllrel1{n} = [cidxs(i) cidxs(j)];
%   end
% end
% 
end

