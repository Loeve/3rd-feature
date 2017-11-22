addpath(genpath(pwd));
% vl_setup

[fname,pname]=uigetfile('*.bmp;*.png;*.jpg','pick picture file','MultiSelect', 'on');
str = fullfile(pname,fname);
image=imread(str);
image = rescale_max_size(image, 200);

% Compute features.
[contrast_map,center_surround_map,color_spatial_map] = ...
    extract_salient_feature_maps(image);
% Infer saliency mask.
 [salient_map,rects] = infer_salient_map(contrast_map, ...
    center_surround_map,color_spatial_map);

showresult(pname,filename{k},rects,image,salient_map, ...
    contrast_map,center_surround_map,color_spatial_map);

 


