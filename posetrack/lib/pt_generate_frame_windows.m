function [windows] = pt_generate_frame_windows(num_frames, win_size, overlap)
% windows: every row contains [startIdx endIdx] of each window

if(nargin < 3)
    overlap = 0;
end

if(num_frames < win_size)
    windows = [1, num_frames];
    return;
end

n = 1;
stIdx = 1;
while (stIdx <= num_frames)
   stIdx = max(1,stIdx-overlap);
   endIdx = min(stIdx+win_size-1, num_frames);
   windows(n, :) = [stIdx, endIdx];
   stIdx = endIdx+1;
   n = n + 1;
end
    
end

