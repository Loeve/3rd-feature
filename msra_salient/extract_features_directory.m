% Extract features from images located in directory
%
% Vicente Ordonez @ Stony Brook University
%                   State University of New York

function [naffected] = extract_features_directory(directory, outputdirectory)
images_list = dir(directory);

disp('Reading and extracting features from images...')
affected = 0;
for i = 1: size(images_list)
    % Check if the item represents an image.
    imagename = [outputdirectory, '/', images_list(i).name];
    if ((images_list(i).isdir == 0) && ...
        isempty(findstr(images_list(i).name, '.ms.jpg')) && ...
        isempty(findstr(images_list(i).name, '.csh.jpg')) && ...
        isempty(findstr(images_list(i).name, '.cwcsd.jpg')) && ...
        (~exist([imagename, '.ms.jpg'], 'file') || ...
         ~exist([imagename, '.csh.jpg'], 'file') || ...
         ~exist([imagename, '.cwcsd.jpg'], 'file')))
    
        imagename = [directory, '/', images_list(i).name];
        
        tic();
        image = imread(imagename);
        [height width channels] = size(image);
        if (width > 400 || height > 400)
            if (width > height)
                height = 400 * height / width;
                width = 400;
            else
                width = 400 * width / height;
                height = 400;
            end
            image = imresize(image, [height width]);
        end
        
        imagename = [outputdirectory, '/', images_list(i).name];

        if (channels == 1)
            image2 = image;
            [height width c] = size(image);
            image = zeros(height, width, 3);
            image(:,:,1) = image2;
            image(:,:,2) = image2;
            image(:,:,3) = image2;
        end
            
        fprintf('Processing for image %s started\n', imagename);
            mscm = multiscale_contrast_map(image);
            imwrite(mscm, [imagename, '.ms.jpg']);

            [map rectangles] = calculate_distances_map(image, 0);
            cshm = center_surround_histogram_map(map, rectangles, 0);
            imwrite(cshm, [imagename, '.csh.jpg']);
            
            [colormaps] = generate_colormaps(image, 6,  0);
            cwcsdm = center_weighted_color_spatial_distribution_map(colormaps);
            imwrite(cwcsdm, [imagename, '.cwcsd.jpg']);
            
            fprintf('Processing image %s took %.2f secs...\n', imagename, toc());
            affected = affected + 1;
    end
end

naffected = affected;

end

