% Takes the input image and calculates the multiscale contrast map of the
% image as defined in 'Learning to Detect a Salient Object' by Liu, Sun,
% Tang & Sum CVPR 2007.
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [multiscale_contrast_map]=multiscale_contrast_map(input_image)

% The number of levels in the scale pyramid.
total_number_of_levels = 6;
[height width channels] = size(input_image);
multiscale_contrast_map = zeros(height, width);

filtered_image = zeros(height, width, channels);

for level = 1 : total_number_of_levels
    % Apply the gaussian reduction to this level of the pyramid.
    for channel = 1 : channels
        filtered_image(:, :, channel) = ...
            apply_gaussian_filter(input_image(:, :, channel));
    end
    
    % Subsample this level of the pyramid
    scaled_image = imresize(filtered_image, 1 / 2 ^ (level - 1), 'bilinear');
    
    % Calculate the contrast map.
    scaled_contrast_map = contrast_map(scaled_image);
    
    % Rescale the contrast map to the original image size.
    rescaled_contrast_map = imresize(scaled_contrast_map, [height width]);
%     figure; imshow(rescaled_contrast_map);
    % imwrite(rescaled_contrast_map, ['testdata/cato' int2str(level) '.jpg']);
    
    % Reassign the image with the current filtered_image for the next level
    % of the pyramid.
    input_image = filtered_image;
    
    multiscale_contrast_map = ...
        multiscale_contrast_map + rescaled_contrast_map;
end

multiscale_contrast_map = ...
    multiscale_contrast_map + max(min(multiscale_contrast_map(:)), 0);
multiscale_contrast_map = ...
    multiscale_contrast_map ./ max(multiscale_contrast_map(:));

% figure; imshow(multiscale_contrast_map);
end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         