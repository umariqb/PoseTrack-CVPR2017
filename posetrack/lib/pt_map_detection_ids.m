function     [pwProbGlobalIds, numDetections] = pt_map_detection_ids(pwProb, detPerFrame);

pwProbGlobalIds = nan(size(pwProb,1), 4);

detStartIndexes = [0; cumsum(detPerFrame(1:end-1))];
numDetections = sum(detPerFrame);

for i = 1:size(pwProb, 1)
    f0 = pwProb(i, 1)+1;
    f1 = pwProb(i, 2)+1;
    d0 = pwProb(i, 3);
    d1 = pwProb(i, 4);
     
    newd0 = detStartIndexes(f0) + d0;
    newd1 = detStartIndexes(f1) + d1;
    
    pwProbGlobalIds(i, 1) = newd0;
    pwProbGlobalIds(i, 2) = newd1;
    pwProbGlobalIds(i, 3) = pwProb(i,5);
    pwProbGlobalIds(i, 4) = pwProb(i,6);
   
end


end

