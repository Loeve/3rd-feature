% Calculate a surrounding rectangle for the given inner rectangle such
% that the area enclosed between the bigger and smaller rectangle is the
% same as the inner rectangle.
% The input inner rectangle is in the format [left top width height]
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [surrounding_rectangle] = ...
    calculate_surrounding_rectangle(inner_rectangle)

left = inner_rectangle(1); top = inner_rectangle(2);
width = inner_rectangle(3); height = inner_rectangle(4);

surrounding_width = floor(sqrt(2) * width);
surrounding_height = floor(sqrt(2) * height);

surrounding_left = left - floor((sqrt(2) - 1) * width / 2);
surrounding_top = top - floor((sqrt(2) - 1) * height / 2);

surrounding_rectangle = ...
   [surrounding_left surrounding_top surrounding_width surrounding_height];
end

