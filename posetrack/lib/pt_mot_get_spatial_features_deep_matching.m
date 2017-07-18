function [X_pos, X_neg] = mot_get_spatial_features_deep_matching(expidx)


RandStream.setGlobalStream ...
        (RandStream('mt19937ar','seed',42));

p = mot_exp_params(expidx);

bVis = false;

save_file = [p.pairwiseDir '/feat_spatial.mat'];

if exist(save_file, 'file') == 2
    fprintf('Loading feature file: %s\n', save_file);
    load(save_file);
    fprintf('Loaded %s\n',save_file);
    return;
end

fprintf('Loading annotations %s\n', p.trainGT);

trainSeq = p.trainSequenceNames;
numSeq = numel(trainSeq);
gtFile = p.trainGT;

annotations = cell(1, numSeq);
for s=1:numSeq
    filePath = sprintf(gtFile, trainSeq{s});
    ann = mot_convert_txt_to_struct(filePath);
    
    if(~isfield(ann,'Xi') || ~isfield(ann,'Yi') || ~isfield(ann,'W') || ~isfield(ann,'H'))
        fprintf('Invalid file %s\n', filePath);
        X_pos = [];
        X_neg = [];
        return;
    end
    annotations{s} = ann;
end


numFrames = 0;
numPersons = 0;
for s=1:numSeq
    ann = annotations{s};
    numFrames  = numFrames  + size(ann.Xi, 1);
    numPersons = numPersons + size(ann.Xi, 2);        
end

dim = 5;
% initialize the arrays, otherwise too slow when using cat
X_pos = zeros(numFrames*numPersons*2*p.max_time_dist,dim,'single');
X_neg = zeros(numFrames*numPersons*2*p.max_time_dist,dim,'single');


% example counter for each class
n_pos = 0;
n_neg = 0;


for s=1:numSeq
    ann = annotations{s};
    numFrames  = size(ann.Xi, 1);
    numPersons = size(ann.Xi, 2);        
    
    for f=1:numFrames
        
        imgPath = fullfile(sprintf(p.imageName, trainSeq{s}, f));
        I1 = imread(imgPath);

        stIdx = max(1,f-p.max_time_dist);
        enIdx = min(numFrames, f+p.max_time_dist);

        for i=stIdx:enIdx
            if(f == i)
                continue;
            end
            
            imgPath = fullfile(sprintf(p.imageName, trainSeq{s}, i));
            I2  = imread(imgPath);
            
            for p = 1:numPersons
                x = ann.Xi(f,p);
                y = ann.Yi(f,p);
                w = ann.W(f,p);
                h = ann.H(f,p);
                bbox1   = [x,y,w,h];
            
                x = ann.Xi(i,p);
                y = ann.Yi(i,p);
                w = ann.W(i,p);
                h = ann.H(i,p);
                bbox2 = [x,y,w,h];

            end
        end
    end

    if (size(X_neg,1) > size(X_pos,1))
        nFeat = size(X_pos,1);
        idxs_rnd = randperm(size(X_neg,1));
        idxsSamp = idxs_rnd(1:nFeat);
        X_neg = X_neg(idxsSamp,:);
    end
    
    [X_norm, opts.X_min, opts.X_max] = getFeatNorm([X_pos;X_neg]);
    X_pos_norm = X_norm(1:size(X_pos,1),:);
    X_neg_norm = X_norm(size(X_pos,1)+1:end,:);
    
    reg_type = 0; % L2
    C = 1e-3;
    
    nPos = size(X_pos_norm,1);
    nNeg = size(X_neg_norm,1);
    lab_pos = ones(nPos,1);
    lab_neg = zeros(nNeg,1);
    lab = [lab_pos; lab_neg];
    ex = sparse(double([X_pos_norm; X_neg_norm]));
    model = train(lab, ex, ['-s ' num2str(reg_type) ' -B 1 -c ' num2str(C)]);

    spatial_model.training_opts = opts;
    spatial_model.log_reg = model;
    save(modelFname, 'spatial_model');
end
% ------------------------------------------------------------------------