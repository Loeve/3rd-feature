function [output_image] = apply_gaussian_filter(input_image)
% Apply a gaussian filter to the input image.

% Create the gaussian convolution kernel and allocate the output_image.
gaussian_kernel = [0.05 0.25 0.4 0.25 0.05];
[height width] = size(input_image);
padding_image = zeros(height + 4, width + 4);

% Copy the value of the input image on the output image.
padding_image(3 : height + 2, 3 : width + 2) = input_image;
padding_image(1, :) = padding_image(3, :);
padding_image(2, :) = padding_image(3, :);
padding_image(height + 3, :) = padding_image(height + 2, :);
padding_image(height + 4, :) = padding_image(height + 2, :);

padding_image(:, 1) = padding_image(:, 3);
padding_image(:, 2) = padding_image(:, 3);
padding_image(:, width + 3) = padding_image(:, width + 2);
padding_image(:, width + 4) = padding_image(:, width + 2);


% Apply convolution in the horizontal direction.
padding_image = conv2(padding_image, gaussian_kernel, 'same');
padding_image = padding_image'; % Prepare for vertical convolution.

% Apply convolution in the vertical direction.
padding_image = conv2(padding_image, gaussian_kernel, 'same');
padding_image = padding_image'; % Restore to the original orientation.

output_image = padding_image(3 : height + 2, 3 : width + 2);

end

