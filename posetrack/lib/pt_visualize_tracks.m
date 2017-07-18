function pt_visualize_tracks( p, detections, labels, vid_dir, num_frames )
%PT_VISUALIZE_TRACKS Summary of this function goes here
%   Detailed explanation goes here

colors = {'r','g','b','c','m','y'};
figure(109); clf;
fn          = dir([vid_dir,'/*.jpg']);

stIdx  = min(detections.frameIndex);
endIdx = max(detections.frameIndex);
for fIdx = stIdx:endIdx
    fprintf('Frame: %d/%d\n',fIdx, num_frames);
    fr_fn = fullfile(vid_dir, fn(fIdx).name);
    frame = imread(fr_fn);
    
    cfIdx = detections.frameIndex == fIdx;
    dets = Detections.slice(detections, cfIdx);
    cLabels = labels(cfIdx);
    
    imagesc(frame); axis equal; hold on;
    
    cb = dets.unPos;
    scale = dets.scale;
    plot(cb(:,1),cb(:,2),'b+');
    for j = 1:size(dets.unPos, 1)
        label = cLabels(j);
        c = colors(mod(label, length(colors))+1);
        text(cb(j,1),cb(j,2), num2str(label), 'Color', c{1}, 'FontSize', 10);
%         text(cb(j,1),cb(j,2)+15, num2str(dets.index(j)), 'Color', c{1}, 'FontSize', 10);
        r = (p.patchSize * 1/scale(j))/2;
        rectangle('Position', [cb(j,1)-r,cb(j,2)-r, 2*r, 2*r], 'EdgeColor', c{1});
    end
    pause();
    hold off;

end

end

