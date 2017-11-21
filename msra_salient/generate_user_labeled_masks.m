% Takes as input the file containing the coordinates of the salient
% rectangles labeled by 9 different users for each image in the
% Microsoft Research dataset for salient objects.
%
% Vicente Ordonez @ Stony Brook University
%                   Staet University of New York

function [success] = generate_user_labeled_masks(filepath)

file_id = fopen(filepath);

% Read the number of items described in the file.
count = textscan(file_id, '%d', 1);
count = count{1};

for i = 1 : count
    % Retrieve the path of the file of the current image.
    image_path = textscan(file_id, '%s', 1);
    
    % Works this way, don't ask why.
    image_path = image_path{1, 1};
    image_path = image_path{1, 1};
    
    % Retrieve the dimensions of the image.
    dimensions = textscan(file_id, '%d %d', 1);
    image_width = dimensions{1, 1};
    image_height = dimensions{1, 2};
    
    % Now retrieve the coordinates of the 9 user labeled rectangles.
    rectangles = textscan(file_id, '%d %d %d %d;', 9);
    lefts = rectangles{1, 1};
    tops = rectangles{1, 2};
    rights = rectangles{1, 3};
    bottoms = rectangles{1, 4};

    % For some reason the masks in the user labeled dataset 
    % go from 0 to image_width including both.
    lefts = 1 + lefts;
    tops = 1 + tops;
    labeled_mask = zeros(image_height, image_width);
    for j = 1: 9
        if (rights(j) >= image_width)
            rights(j) = rights(j) - 1;
        end
        
        if (bottoms(j) >= image_height)
            bottoms(j) = bottoms(j) - 1;
        end
        labeled_mask(tops(j) : bottoms(j), lefts(j) : rights(j)) = ...
            labeled_mask(tops(j) : bottoms(j), lefts(j) : rights(j)) + 1;
    end
    labeled_mask = labeled_mask ./ 9;
 
    imwrite(labeled_mask, ['Image\' image_path '.label.jpg']);
    
    fprintf('Saving mask %d...\n', i);
end

success = rectangles;
end

