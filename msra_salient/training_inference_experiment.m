function training()

% FIRST LOAD ALL THE USER LABELED MASKS INTO THE APPROPRIATE STRUCTURE.
% Load Training Data.
directory_path = 'Image/0';
output_directory = 'recognition4';
trainNdx = 1 : 10;
testNdx = 1 : 5;

% List all the images in the training images directory.
images_list = dir(directory_path);

scale = 50;
training_set_length = size(trainNdx, 2);
test_set_length = size(testNdx, 2);
labeled_masks =  zeros(scale, scale, training_set_length);
test_images = zeros(scale, scale, test_set_length);
labeled_masks_test = zeros(scale, scale, test_set_length);
feature_maps = zeros(1, scale, scale, training_set_length);
feature_maps_test = zeros(1, scale, scale, test_set_length);

% Load the training labeled masks.
% Load the feature-map features for some
% traning and test images.
count = 0;
iterator = 1;
while count < training_set_length + test_set_length
    if ((images_list(iterator).isdir == 0) && ...
         ~isempty(findstr(images_list(iterator).name, '.label.jpg')))
        count = count + 1;
        image_path = [directory_path, '/', images_list(iterator).name];
        mask = imread(image_path);
        
        % Calculate the threshold to produce a binary image.
        % threshold = max(mask(:));%graythresh(mask);
        
        % Make the image binary.
        training_mask = binarize_mask(mask);
        
        % training_mask = binarize_labeled_mask(mask);
        % training_mask = 2 * training_mask - 1;
        % Now make the binary image [1, -1].
        training_mask = imresize(training_mask, [scale scale]);
        training_mask = sign(double(training_mask) - 1);
        %training_mask
        
        % Store the labeled mask.
        if (count <= training_set_length)
            labeled_masks(:, :, count) = training_mask;
        else 
            labeled_masks_test(:, :, count - training_set_length) = training_mask;
            %figure;imshow(labeled_masks_test(:, :, count - training_set_length));
        end
        
        % Calculate the name of the corresponding file with the
        % pre-calculated feature map.
        feature_map_path = ...
            regexprep(images_list(iterator).name, '.label.jpg', '.csh.jpg');
        
        image_path = ...
            regexprep(images_list(iterator).name, '.label.jpg', '');
        image_path = strcat(directory_path, '/', image_path);
        
        if (count > training_set_length)
            image_paths_test{count - training_set_length} = image_path;
        end
        
        feature_map_path = [directory_path, '/', feature_map_path];
        feature_map = imread(feature_map_path);
        feature_map = imresize(feature_map, [scale scale]);
        feature_map = im2double(feature_map);
        feature_map = 8.0 * feature_map - 3.0; % For center surround hist.
        %feature_map = 10.0 * feature_map - 2.0; % For multiscale contrast.
        %feature_map = 8.0 * feature_map - 3.0; % For cs color histograms.
        
        if (count <= training_set_length)
            feature_maps(1, :, :, count) = feature_map;
        else
            image = imread(image_path);
            image = imresize(image, [scale scale]);
            image = rgb2gray(image);
            feature_maps_test(1, :, :,count - training_set_length) = ...
                feature_map;
            test_images(:, :, count - training_set_length) = image;
        end
        
    end
    iterator = iterator + 1;
end


figure;
for i = testNdx
    subplot(2,test_set_length,i);
    lb = permute(feature_maps_test(1, :,:,i), [2 3 4 1]);
    imshow(lb);
end

for i = test_set_length + 1 : 2 * test_set_length
    subplot(2,test_set_length,i);
    imshow(labeled_masks_test(:,:,i - test_set_length));
    labeled_masks_test(:,:,i - test_set_length);
end
suptitle('test data')

% I will figure it out later what nstates mean.
nstates = 2;

% Make features and feature engine.
featureEng = latticeFeatures(0, 0);

% Supply the feature maps to the CRF framework.
trainFeatures = feature_maps;
traindata.nodeFeatures = mkNodeFeatures(featureEng, trainFeatures);
traindata.edgeFeatures = mkEdgeFeatures(featureEng, trainFeatures);
traindata.nodeLabels = labeled_masks;
traindata.ncases = length(trainNdx);
trainNdx = 1:traindata.ncases;

nNodeFeatures = size(traindata.nodeFeatures, 1);
nEdgeFeatures = size(traindata.edgeFeatures, 1);
winit = initWeights(featureEng, nNodeFeatures, nEdgeFeatures);

testFeatures = feature_maps_test;
testdata.nodeFeatures = mkNodeFeatures(featureEng,testFeatures);
testdata.edgeFeatures = mkEdgeFeatures(featureEng,testFeatures);
testdata.nodeLabels = labeled_masks_test;
testdata.ncases = length(testNdx);
testNdx = 1:testdata.ncases;

% Random params
infEng = latticeInferBP(nstates);
% showResults(winit, testNdx, featureEng, infEng, testdata, 'rnd params');


% BFGS training with Belief - Propagation
reg = 1;
maxIter = 3;
options = optimset('Display','iter','Diagnostics','off','GradObj','on',...
    'LargeScale','off','MaxFunEval',maxIter);

gradFunc = @scrfGradient;
gradArgs = {featureEng, infEng, traindata, reg};

weights = fminunc(gradFunc,winit,options,trainNdx,gradArgs{:});
weights

showResults(weights, testNdx, featureEng, infEng, testdata, 'BFGS+BP', image_paths_test, output_directory);


end

function showResults(weights, testNdx, featureEng, infEng, testdata, ttl, image_paths_test, output_directory)

%figure;
for i = testNdx
    % subplot(2,5,i);
    featureEng = enterEvidence(featureEng, testdata, i);
    [nodePot, edgePot] = mkPotentials(featureEng, weights);
    [nodeBel, MAPlabels] = infer(infEng, nodePot, edgePot);
    
    original_image = imread(image_paths_test{i});
    ground_truth = imread([image_paths_test{i} '.label.jpg']);
    image_name = regexprep(image_paths_test{i}, 'Image/0/', '');
    image_name
    imwrite(original_image, [output_directory '/' image_name]); 
    imwrite(ground_truth, [output_directory '/' image_name '.ground.jpg']);
    [height width channels] = size(original_image);
    crf_inference_map = imresize(MAPlabels, [height width]);
    imwrite(crf_inference_map, [output_directory '/' image_name '.feature3.jpg']);
    
    
    % imshow(MAPlabels);
    %title(sprintf('%d',i));
end
%testErrorRate = classifPerformance(weights, testNdx, featureEng, infEng, testdata);
%suptitle(sprintf('%s, error rate = %5.3f', ttl, testErrorRate))
%drawnow
end

function [output_mask] =  binarize_mask(mask)
    threshold = (max(mask(:))) / 2;
    output_mask = zeros(size(mask, 1), size(mask, 2));
    for i = 1 : size(mask, 1)
        for j = 1: size(mask, 2)
            if (mask(i, j) >= threshold)
                output_mask(i, j) = 255;
            else
                output_mask(i, j) = -255;
            end
        end
    end
end
