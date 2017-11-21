% Generate a rectangle in the format [left top width height]
% such that this rectangle is centered at (x,y), the smallest size is
% of size scale and the aspect ratio is the given input aspect_ratio.
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [rectangle] = generate_rectangle(x, y, scale, aspect_ratio)

% Calculate the current rectangle size and location.
if (aspect_ratio > 1.0)
    rectangle_height = scale;
    rectangle_width = floor(rectangle_height / aspect_ratio);
else
    rectangle_width = scale;
    rectangle_height = floor(rectangle_width * aspect_ratio);
end

rectangle_left = x - floor(rectangle_width / 2);
rectangle_top = y - floor(rectangle_height / 2);
rectangle = [rectangle_left, rectangle_top,...
             rectangle_width, rectangle_height];
end

