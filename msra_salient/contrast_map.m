% Calculates a contrast map using neighbouring pixels.
% The input image is expected to be a rgb image.
%
% Author: Vicente Ordonez @ Stony Brook University
%                           State University of New York

function [contrast_map] = contrast_map(input_image)

dx = [-1 0 1; -1 0 1; -1 0 1];
dy = [-1 -1 -1; 0 0 0; 1 1 1];

[height,width,channels] = size(input_image);

% Initialize the output buffer image.
contrast_map = zeros(height, width);

% Traverse the image using a sliding window of 9x9.
for row = 2 : height - 1
    for column = 2 : width - 1
        contrast_value = 0;
        for window_index = 1 : 9
           % Vector representing the intensity of the current pixel.
           intensity = input_image(row, column, :); 
           % Vector representing the intensity of the current
           % neighbouring pixel.
           neighbor_intensity = input_image(row + dy(window_index),... 
                                            column + dx(window_index), :);
           
           % diff = intensity - neighbor_intensity
           % Use the euclidean distance as a measure of distance.
           euclidean_distance_squared =...
               sum((intensity - neighbor_intensity) .^ 2);
           
           % Accumulate the distances to all pixels in the neighborhood.
           contrast_value = contrast_value + euclidean_distance_squared;
        end
        contrast_map(row, column) = contrast_value;
    end
end

% Rescale pixel values.
contrast_map = contrast_map + max(min(contrast_map(:)), 0);
contrast_map = contrast_map ./ max(contrast_map(:));

end

