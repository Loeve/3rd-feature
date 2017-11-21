% Calculate the spatial distribution map from the colormaps obtained
% using the function generate_colormaps_using_GMM.
%
% Author: Vicente Ordonez @ Stony Brook University
%                           State University of New York

function [output_map] = color_spatial_distribution_map(colormaps)

[height width numberOfColors] = size(colormaps);
output_map = zeros(height, width);

variance = calculate_spatial_variances(colormaps);

for color = 1: numberOfColors
   output_map = ...
       output_map + (1 - variance(color)) * colormaps(:, :, color);
end

output_map = output_map ./ max(output_map(:));

end