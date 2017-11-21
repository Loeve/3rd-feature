% Calculates the saliency of the pixel (row, column) and the most salient
% rectangle for this position.
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [saliency_value,most_salient_rectangle] = ...
    saliency_at_position(input_image, integral_histogram, row, column, debug)

% The range of the rectangle size will be [0.1, 0.7] * min(width, height).
[height,width,channels] = size(input_image);
%pad_horizontal = min(column, width - column);
%pad_vertical = min(row, height - row);
%pad = min(pad_horizontal, pad_vertical);

% Calculating the range of sizes of the candidate rectangles.
initial_side_length = floor(0.1 * min(height, width)); % floor(0.1 * min(height, width));
final_side_length = floor(0.7 * min(height, width)); %min(sqrt(2) * pad, floor(0.7 * min(height, width)));
length_step = floor(0.1 * min(height, width));

if (debug >= 1)
    imshow(input_image);
    impixelinfo
end

best_rectangle = [0 0 0 0];
% The range of the aspect rations of the rectangles.
aspect_ratios = [0.5 0.75 1.0 1.5 2.0];
max_chi_square_distance = 0;
counter = 0;
    for scale = initial_side_length: length_step : final_side_length
        for aspect_index = 1: length(aspect_ratios)
            aspect_ratio = aspect_ratios(aspect_index);
 
            target_rectangle = ...
                generate_rectangle(column, row, scale, aspect_ratio);
            surrounding = calculate_surrounding_rectangle(target_rectangle);
            
            % Checking that surrounding rectangle falls inside the image.
            if (check_that_rectangle_is_inside(surrounding, width, height)) 
            
                histogram_rectangle = ...
                    calculate_histogram_fast(input_image, ...
                                             integral_histogram, ...
                                             target_rectangle);

                histogram_enclosing_rectangle = ...
                    calculate_histogram_fast(input_image, ...
                                             integral_histogram, ...
                                             surrounding);

                histogram_surrounding = ...
                    histogram_enclosing_rectangle - histogram_rectangle;

                distance = chi_square_distance(histogram_rectangle, ...
                                               histogram_surrounding);

                if (debug >= 2)
                    rectangle('Position', target_rectangle, 'EdgeColor', [0 1 0]);
                    counter = 1 + counter;
                    fprintf('Rectangle(%d) = %f\n', counter, distance);
                end

                if (distance > max_chi_square_distance)
                    max_chi_square_distance = distance;
                    best_rectangle = target_rectangle;
                end
            end
        end
    end

if (final_side_length >= initial_side_length)
    saliency_value = max_chi_square_distance;
    most_salient_rectangle = best_rectangle;
    if (debug >= 1) 
        surrounding = calculate_surrounding_rectangle(most_salient_rectangle);
        rectangle('Position', most_salient_rectangle, 'EdgeColor', [1 0 0], 'LineWidth', 2);
        rectangle('Position', surrounding, 'EdgeColor', [0 0 1]);
    end
else
    most_salient_rectangle = [0 0 0 0];
    saliency_value = 0;
end

    function [inside] = check_that_rectangle_is_inside(rectangle, imwidth, imheight)
        rleft = rectangle(1); rtop = rectangle(2);
        rwidth = rectangle(3); rheight = rectangle(4);
        inside = true;
        if (rleft < 0) inside = false; end
        if (rtop < 0) inside = false; end
        if (rleft + rwidth >= imwidth) inside = false; end
        if (rtop + rheight >= imheight) inside = false; end
    end

end

% Test using these commands
% venado = imread('testdata/tester.jpg')
% saliency_at_position(venado, 185, 240)
% saliency_at_position(venado, 95, 550)

