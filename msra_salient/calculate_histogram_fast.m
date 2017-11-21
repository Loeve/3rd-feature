% Calculates the histogram of the section of the input_image overlapped
% with the given rectangle. The input is assumed to be a single channel
% intensity image.
% The input rectangle is in the format [left top width height].
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [histogram] = calculate_histogram_fast(input_image, ...
                                                integral_image, ...
                                                rectangle)

[image_height,image_width,channels] = size(input_image);
image_rectangle = [1 1 image_width image_height];
% encompass°üÎ§£¬Î§ÈÆ
[left,top,width,height] = encompass(rectangle, image_rectangle);

iMin = top;
iMax = top + height - 1;

jMin = left;
jMax = left + width - 1;

boxes = zeros(1, 4);

boxes(1, :) = [iMin jMin iMax jMax];
boxes = uint32(boxes);


histograms = vl_samplinthist(integral_image, boxes);

histogram = histograms(:, 1);