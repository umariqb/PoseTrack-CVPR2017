function [ unary_maps, locreg_pred, nextreg_pred, rpn_prob, rpn_bbox ] = extract_features_multiscale( im, net, p, annorect, pad_orig, pairwise, crop )

scales = p.scales;

fprintf('Extracting multiscale features.\n');
for i=1:length(scales)
    scale_factor = scales(i);
    fprintf('Scale %d/%d (%f)\n', i, length(scales), scale_factor);
    p.scale_factor = scale_factor;
    [unary_maps{i}, locreg_pred{i}, nextreg_pred{i}, rpn_prob{i}, rpn_bbox{i}] = ...
                        extract_features(im, net, p, [], pad_orig, pairwise, crop);
end

