% This function constructs a labeled_map of rgb values.
% It assigns a label from 10 possible labels to each of the rgb components.
% Then it creates a single label accross all the 3 rgb components and 
% returns this map of labels.
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [labeled_map] = rgb_labeled_map(rgb_image)

% Fix to work with any type of images.
rgb_image = im2uint8(rgb_image);

% Each channel is asiggned to a 0 to 9 value.
red_channel = idivide((uint32(99 * ((double(rgb_image(:, :, 1)) / 255)))), 10);
green_channel = idivide((uint32(99 * ((double(rgb_image(:, :, 2)) / 255)))), 10);
blue_channel = idivide((uint32(99 * ((double(rgb_image(:, :, 3)) / 255)))), 10);

% A labeled map is constructed where each values goes from 1 to 1000.
labeled_map = 100 * red_channel  + 10 * green_channel + blue_channel + 1;

end

