% Wrapper function that extracts the saliency map
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [saliency_map, saliency_map2, energy] = extract_saliency_map(image)
    [cm cshm csdm] = extract_salient_feature_maps(image);
    [mask mask2 value] = infer_salient_map(cm, cshm, csdm);
    saliency_map = mask;
    saliency_map2 = mask2;
    energy = value;
end
