function  pt_train_spatial_part_regr(expidx, cidx)

fprintf('train_spatial_part_regr()\n');

if (ischar(expidx))
    expidx = str2num(expidx);
end

if (ischar(cidx))
    cidx = str2num(cidx);
end

p = pt_exp_params(expidx);

cache_dir = fullfile(p.exp_dir, 'scoremaps');
fprintf('loading scoremaps from: %s\n', cache_dir);

ptPairwiseDir = p.ptPairwiseDir;
if (~exist(ptPairwiseDir,'dir')),mkdir(ptPairwiseDir);end
modelFname = [ptPairwiseDir '/spatial_model_cidx_' num2str(cidx)];

try
    %assert(false);
    load(modelFname,'spatiall_model');
    fprintf('cidx: %d - %d\n',cidxs);
    fprintf('spatial model file loaded. quitting.\n');
catch
    [X_pos, X_neg] = pt_get_spatial_features_dm(expidx,cidx);
    
    opts.X_pos_mean = mean(X_pos);
    opts.X_pos_std = std(X_pos);
    idxs1 = X_pos(:,1) >= opts.X_pos_mean(:,1) - 3*opts.X_pos_std(:,1) & X_pos(:,1) <= opts.X_pos_mean(:,1) + 3*opts.X_pos_std(:,1);
    idxs2 = X_pos(:,2) >= opts.X_pos_mean(:,2) - 3*opts.X_pos_std(:,2) & X_pos(:,2) <= opts.X_pos_mean(:,2) + 3*opts.X_pos_std(:,2);
    X_pos = X_pos(idxs1 & idxs2,:);
    
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
    C = 1;
    epsilon = 0;
    
    
    nPos = size(X_pos_norm,1);
    nNeg = size(X_neg_norm,1);
    lab_pos = ones(nPos,1);
    lab_neg = zeros(nNeg,1);
    lab = [lab_pos; lab_neg];
    ex = sparse(double([X_pos_norm; X_neg_norm]));
    model = train(lab, ex, ['-s ' num2str(reg_type) ' -B 1 -c ' num2str(C) '-e ' num2str(epsilon)]);

    temporal_model.training_opts = opts;
    temporal_model.log_reg = model;
    save(modelFname, 'temporal_model');
    
    if(1)
        [pred_lab,acc,prob_out] = predict(lab, ex, model, '-b 1');         
    end
    
    if (0)
        visDir = [pairwiseDir '/vis/'];
        if (~exist(visDir,'dir')),mkdir(visDir);end
        [~,acc,pred] = predict(lab, ex, model, '-b 1');
        scrsz = get(0,'ScreenSize');
        figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
        vis_logreg(pred,acc,1:nPos,1+nPos:nPos+nNeg);
        print(gcf,'-dpng',[visDir '/logreg_cidx_' num2str(cidxs(1)) '_' num2str(cidxs(2)) '.png']);
        close all;
    end
end


end

