% Calculate the spatial distribution map from the colormaps obtained
% using the function generate_colormaps_using_GMM. This additionally
% weights less those color areas close to the image boundaries.
% 
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [output_map] = ...
    center_weighted_color_spatial_distribution_map(colormaps)


[height,width,numberOfColors] = size(colormaps);
output_map = zeros(height, width);

variance = calculate_spatial_variances(colormaps);

% Calculate the weights D(c) described in the paper 'Learning to Detect a 
% Salient Object' by Liu et al. CVPR 2007.
center_distance_weights = zeros(numberOfColors, 1);
pixels_in_cluster = 0;
epsilon = 0.00001;
for color = 1: numberOfColors
    for y = 1: height
        for x = 1: width
            color_probability = colormaps(y, x, color);
            if (color_probability > 0 + epsilon)
                pixels_in_cluster = pixels_in_cluster + 1;
            end
            
            deltay = y - height / 2;
            deltax = x - width / 2;
            distance = sqrt(deltay * deltay + deltax * deltax);
            center_distance_weights(color) = ...
                center_distance_weights(color) + ...
                    colormaps(y, x, color) * distance;
        end
    end
    % Normalize the distance weight. (This is not in the original paper).
    %center_distance_weights(color) = ...
     %   center_distance_weights(color) / pixels_in_cluster;
    %pixels_in_cluster = 0;
end

% Normalize in the range [0-1].
center_distance_weights = ...
    (center_distance_weights - min(center_distance_weights(:))) ./ ...
    (max(center_distance_weights(:)) - min(center_distance_weights(:)));

for color = 1: numberOfColors
   output_map = ...
       output_map + (1 - variance(color)) * ...
        (1 - center_distance_weights(color)) * colormaps(:, :, color);
end

output_map = output_map ./ max(output_map(:));

end