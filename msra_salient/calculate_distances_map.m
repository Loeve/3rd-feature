function [map,rectangles] = calculate_distances_map(image, debug)
[height,width,~] = size(image);
map = zeros(height, width);
rectangles = zeros(height, width, 4);

% The labeled map has each pixel labeled from 1 to 1000.
labeled_map = rgb_labeled_map(image);

% Calculate the integral histogram for fast calculations.
integral_histogram = vl_inthist(uint32(labeled_map));

step = floor(max(5, max(height, width) / 36));

% Store the minimum distance found in the minimum variable.
minimum = flintmax;  % Initialize with the maximum double value.
for i = 1 : step : height
    for j = 1 : step : width
        [distance,rectangle] = ...
            saliency_at_position(labeled_map, integral_histogram, i, j, 0);
        map(i, j) = distance;
        if (distance ~= 0 && distance < minimum)
            minimum = distance;
        end
        rectangles(i, j, :) = rectangle;
        if (debug >= 1)
            fprintf('(%d, %d) = %d\n', i, j, map(i,j));
        end
    end
end

% Normalize the distances map to the [0 ... 1] range.
for i = 1: step : height
    for j = 1: step : width
        if (map(i,j) ~= 0)
            map(i,j) = map(i,j) - minimum;
        end
    end
end
map = map / max(map(:));