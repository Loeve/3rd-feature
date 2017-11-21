% Use a CRF to infer a saliency mask based on the input feature maps.
% Based on the method sketched on Liu et al. 'Learning to Detect a Salient Object'
% and using the K. Murphy CRF for 2D Lattices library.
%
% Author: Vicente Ordonez @ Stony Brook University
%                           State University of New York

function [salient_map,rect] = infer_salient_map(contrast_map, center_surround_map, color_spatial_map)

[height,width,c] = size(contrast_map);
% scale = 100;   取消了尺寸变换

% contrast_map = im2double(imresize(contrast_map, [scale scale]));
contrast_map = 10.0 * contrast_map - 2.0; % For multiscale contrast.

% center_surround_map = im2double(imresize(center_surround_map, [scale scale]));
center_surround_map = 8.0 * center_surround_map - 3.0; % For center surround hist.

% color_spatial_map = im2double(imresize(color_spatial_map, [scale scale]));
color_spatial_map = 8.0 * color_spatial_map - 3.0; % For cs color histograms.

% Prepare feature maps for inference into CRF framework.
feature_maps = zeros(3, height,width, 1);
feature_maps(1, :, :, 1) = contrast_map;
feature_maps(2, :, :, 1) = center_surround_map;
feature_maps(3, :, :, 1) = color_spatial_map;

% Prepare feature engine.
featureEng = latticeFeatures(0, 0);

nstates = 2;
testFeatures = feature_maps;
testdata.nodeFeatures = mkNodeFeatures(featureEng,testFeatures);
testdata.edgeFeatures = mkEdgeFeatures(featureEng,testFeatures);
labeled_masks_test = zeros(height,width, 1); % Dummy node labels.
testdata.nodeLabels = labeled_masks_test;
testdata.ncases = 1;

% This weight parameters were learned using 2000 images and user labeled
% masks from the MSRA Salient Object Database.
weights =[-0.1245 0.5379 0.7741 0.7778 0.8486 0.2229 0.3007 0.3843]; 

featureEng = enterEvidence(featureEng, testdata, 1);
[nodePot, edgePot] = mkPotentials(featureEng, weights);

infEng = latticeInferBP(nstates);
[nodeBel, MAPlabels, niter, edgeBel, logZ] = infer(infEng, nodePot, edgePot);   %MAPlabels：-1,1二值图


crf_inference_map=(MAPlabels+1)*0.5;     %  0,1二值图

salient_map = crf_inference_map;
% output_energy = logZ;

bw_img=logical(salient_map);
% T = graythresh(log_sal_map);  
% bw_img = im2bw(salient_map, T);
img_reg= regionprops(bw_img,  'area', 'boundingbox');  
areas= [img_reg.Area];  
rects= cat(1,  img_reg.BoundingBox);  
[~, max_id] = max(areas);  
rect = rects(max_id, :);  
 

end
