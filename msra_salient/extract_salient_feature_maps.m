% Function that extracts the feature maps described in the CVPR 2007
% paper 'Learning to Detect a Salient Object' by Microsoft Research Asia.
%
% Author: Vicente Ordonez @ Stony Brook University
%                           State University of New York

function [contrast_map, center_surround_map, color_spatial_map] = ...
    extract_salient_feature_maps(input_image)
    
    input_image = rescale_max_size(input_image, 400);
    %if black and white
    if(size(input_image,3) < 3)
        %black and white image
        input_image = cat(3, input_image, input_image, input_image); %make it a trivial color image
    end

    % Calculating the multiple resolution contrast map. 
    contrast_map = multiscale_contrast_map(input_image);
    
    % Calculating the Center surround histogram map.
    [map,rectangles] = calculate_distances_map(input_image, 0);
    center_surround_map = center_surround_histogram_map(map, rectangles, 0);
    
    % Calculating the Center weighted color spatial distribution map.
    [colormaps] = generate_colormaps(input_image, 6,  0);
    color_spatial_map = center_weighted_color_spatial_distribution_map(colormaps);

end
