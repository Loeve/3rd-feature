% Center Surround Histogram Feature from the paper 'Learning to Detect a Salient Object'
% by Liu et. al CVPR 2007. 
%
% Author: Vicente Ordonez @ Stony Brook University
%                           State University of New York

function [map] = center_surround_histogram_map(distances, rectangles, debug)
% At this point we have calculated the most distinctive rectangle for
% each pixel in the input image. Now proceed to calculate the output_map.
[height,width] = size(distances);
map = zeros(height, width);

step = floor(max(5, max(height, width) / 36));

for i = 1: step: height 
    for j = 1: step: width
        
        % Pick the most distinctive rectangle for the current position.
        most_distinctive_rectangle = rectangles(i, j, :);
        
        % Fix the rectangle so that it fits inside the image limits.
        [r_left,r_top,r_width,r_height] = ...
           encompass(most_distinctive_rectangle, [1 1 width height]);
        
        % Get the chi square distance for this pixel.
        chi_square_at_this_pixel = distances(i, j);
        
        % Inverse of the variance for the gaussian fall off.
        inverse_sigma_squared = 3 / (r_height * r_width);
        
        for idy = r_top: r_top + r_height - 1
            for idx = r_left: r_left + r_width - 1
                
              % Calculate the weight fall off for the current element.
              deltay = idy - i; deltax = idx - j;
              euclidean_distance_squared = deltay * deltay + deltax * deltax;
              
              %sigma_squared = max(r_height, r_width) / 3
              weight = exp(-0.5 * euclidean_distance_squared * inverse_sigma_squared);
              % weight = 1;
              % chi_square_at_this_pixel = 1;
              % Acumulate the weighted chi square distances on the output_map.
              map(idy, idx) = ...
                  map(idy, idx) + weight * chi_square_at_this_pixel;
            end
        end
        
    end
end

map = map / max(map(:));

if (debug > 0)
    imshow(map)
    impixelinfo
end

end

