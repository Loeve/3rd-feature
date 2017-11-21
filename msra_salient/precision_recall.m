% Calculates the precision and recall for the test data located in
% the direction specified by the input directory parameter.
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [precision, recall, F] = precision_recall(directory)

images_list = dir(directory);


count = 0;
total_detected = 0;
total_ground_truth = 0;
correctly_detected = 0;
for i = 1: size(images_list)
    if ((images_list(i).isdir == 0) && ...
        ~isempty(findstr(images_list(i).name, '.ground.jpg')))
        count = count + 1;
        ground_truth_path = [directory '/' images_list(i).name];
        inferred_map_path = ...
            regexprep(images_list(i).name, '.ground.jpg', '.infer.jpg');
        inferred_map_path = [directory '/' inferred_map_path];
        
        ground_truth = im2double(imread(ground_truth_path));
        inferred_map = im2double(imread(inferred_map_path));
        
        [left top right bottom] = enclosing_rectangle(inferred_map);
        
        enclosing_map = zeros(size(inferred_map, 1), size(inferred_map, 2));
        if (top < bottom && left < right)
            enclosing_map(top:bottom, left:right) = ...
                1 + enclosing_map(top:bottom, left:right);
        else
            ground_truth_path
        end
        
        [height width] = size(enclosing_map);
        ground_truth = imresize(ground_truth, [height width], 'nearest');
        if (size(ground_truth, 1) == size(enclosing_map, 1) && ...
            size(ground_truth, 2) == size(enclosing_map, 2))
        
            correctly_detected = ...
                correctly_detected + sum(sum(ground_truth .* enclosing_map)); 
        else
            ground_truth_path
            size(ground_truth)
            size(enclosing_map)
        end
        
        total_detected = total_detected + sum(enclosing_map(:));
        total_ground_truth = total_ground_truth + sum(ground_truth(:));
        
        %figure;imshow([inferred_map enclosing_map]);
        
    end
end

precision = correctly_detected / total_detected;
recall = correctly_detected / total_ground_truth;
F = (1.5 * precision * recall) / (0.5 * precision + recall);

end


function [left top right bottom] = enclosing_rectangle(map)
    [height width] = size(map);
    top = -1;
    bottom = -1;
    left = -1;
    right = -1;
    for i = 1: height
        for j = 1: width
            if (map(i, j) >= 1 && top < 0)
                top = i;
            end
            if (map(i,j) >= 1)
                bottom = i;
            end
        end
    end
    
    for i = 1: width
        for j = 1: height
            if (map(j, i) >= 1 && left < 0)
                left = i;
            end
            if (map(j, i) >= 1)
                right = i;
            end
        end
    end
    
end
