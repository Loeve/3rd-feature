addpath(genpath('/home/vicente/saliency/CRF2D'));

% FIRST LOAD ALL THE USER LABELED MASKS INTO THE APPROPRIATE STRUCTURE.
% Load Training Data.
%directory_path = 'C:\aware\Datasets\Photonet';
%output_directory = 'C:\aware\Datasets\Photonet_output2';
directory_path = ['/mnt/vicente/data_dp_challenge/saliency/feature_maps/' prefix];

% No need for output directory.
output_directory = ['/mnt/vicente/data_dp_challenge/saliency/saliency_masks/' prefix];

trainNdx = 1 : 1;

% 800 images per folder.
testNdx = 1 : 800;

% List all the images in the images directory.
images_list = dir(['/mnt/vicente/data_dp_challenge/' prefix]);

scale = 50;
training_set_length = 0; %size(trainNdx, 2);
test_set_length = size(testNdx, 2);

test_images = zeros(scale, scale, test_set_length);

labeled_masks =  zeros(scale, scale, training_set_length);
labeled_masks_test = zeros(scale, scale, test_set_length);

feature_maps = zeros(3, scale, scale, training_set_length);
feature_maps_test = zeros(3, scale, scale, test_set_length);

% Load the training labeled masks.
% Load the feature-map features for some
% traning and test images.
count = 0;
iterator = 1;

fprintf('Loading images...');
tic();
while count < training_set_length + test_set_length
    if ((images_list(iterator).isdir == 0) && ...
         isempty(findstr(images_list(iterator).name, '.jpg.cwcsd.jpg')) && ...
         isempty(findstr(images_list(iterator).name, '.jpg.ms.jpg')) && ...
         isempty(findstr(images_list(iterator).name, '.jpg.csh.jpg')))
        fprintf('.');
        if (mod((int32(count)), 50) == 0)
            fprintf('\n');
        end
        count = count + 1;
        
        image_path = [directory_path, '/', images_list(iterator).name];
        image_path = images_list(iterator).name;
        image_path = strcat(directory_path, '/', image_path);

        % Calculate the name of the corresponding file with the
        % pre-calculated feature map. ms = multiscale contrast.
        feature_map_path_ms = ...
            regexprep(images_list(iterator).name, '.jpg', '.jpg.ms.jpg');
        feature_map_path_ms = [directory_path, '/', feature_map_path_ms];
        
        if (exist(feature_map_path_ms, 'file'))
            
            if (count > training_set_length)
                image_paths_test{count - training_set_length} = ...
                    ['/mnt/vicente/data_dp_challenge/' prefix '/', images_list(iterator).name];
            end
            
            feature_map_ms = imread(feature_map_path_ms);
            feature_map_ms = imresize(feature_map_ms, [scale scale]);
            feature_map_ms = im2double(feature_map_ms);

            % Calculate the name of the corresponding file with the
            % pre-calculated feature map. csh = center surround histogram.
            feature_map_path_csh = ...
                regexprep(images_list(iterator).name, '.jpg', '.jpg.csh.jpg');

            feature_map_path_csh = [directory_path, '/', feature_map_path_csh];
            feature_map_csh = imread(feature_map_path_csh);
            feature_map_csh = imresize(feature_map_csh, [scale scale]);
            feature_map_csh = im2double(feature_map_csh);

            % Calculate the name of the corresponding file with the
            % pre-calculated feature map. cwcsd = center weighted color space.
            feature_map_path_cwcsd = ...
                regexprep(images_list(iterator).name, '.jpg', '.jpg.cwcsd.jpg');

            feature_map_path_cwcsd = [directory_path, '/', feature_map_path_cwcsd];
            feature_map_cwcsd = imread(feature_map_path_cwcsd);
            feature_map_cwcsd = imresize(feature_map_cwcsd, [scale scale]);
            feature_map_cwcsd = im2double(feature_map_cwcsd);

            % Scale the values of the features to conform to the CRF library.
            feature_map_ms = 10.0 * feature_map_ms - 2.0; % For multiscale contrast.
            feature_map_csh = 8.0 * feature_map_csh - 3.0; % For center surround hist.
            feature_map_cwcsd = 8.0 * feature_map_cwcsd - 3.0; % For cs color histograms.

            if (count <= training_set_length)
                feature_maps(1, :, :, count) = feature_map_ms;
                feature_maps(2, :, :, count) = feature_map_csh;
                feature_maps(3, :, :, count) = feature_map_cwcsd;
            else
                index = count - training_set_length;
                feature_maps_test(1, :, :, index) = feature_map_ms;
                feature_maps_test(2, :, :, index) = feature_map_csh;
                feature_maps_test(3, :, :, index) = feature_map_cwcsd;
            end
        else
           count = count - 1;
        end
        
    end
    iterator = iterator + 1;
end
endtime = toc()
fprintf('\n\n All images, labels  and feature maps loaded in %d.\n\n', endtime);
fprintf('\n\nTotal images %d\n\n', count);


fprintf('Accomodating the features and labels into the appropriate data structures...\n');
% I will figure it out later what nstates mean.
nstates = 2;

% Make features and feature engine.
featureEng = latticeFeatures(0, 0);

fprintf('Constructing the CRF...\n');
tic();
testFeatures = feature_maps_test;
testdata.nodeFeatures = mkNodeFeatures(featureEng,testFeatures);
testdata.edgeFeatures = mkEdgeFeatures(featureEng,testFeatures);
testdata.nodeLabels = labeled_masks_test;
testdata.ncases = length(testNdx);
testNdx = 1:testdata.ncases;
endtime = toc();
fprintf('CRF construction took %d.\n', endtime);

fprintf('Done with preparing data.\n');

% Random params
infEng = latticeInferBP(nstates);

fprintf('Now let start training using Belief Propagation...\n');
tic();

% Precomputed weights obtained during training.
%weights =[-0.1366 0.3905 0.7304 0.7681 0.1053 0.1664 0.1232 0.3112]; 1000
weights =[-0.1245 0.5379 0.7741 0.7778 0.8486 0.2229 0.3007 0.3843]; %2000
    
weights

fprintf('\n Now let start doing inference and saving results...');
%showResults(weights, testNdx, featureEng, infEng, testdata, 'BFGS+BP', image_paths_test, directory_path, output_directory);

%fprintf('\n Training and detection all completed successfully...');

%function showResults(weights, testNdx, featureEng, infEng, testdata, ttl, image_paths_test, input_directory, output_directory)

ttl = 'BFGS+BP';
input_directory = directory_path;



for i = testNdx
    
    %if (exist(already_processed{i}, 'file'))
    %    continue;
    %end
    
    featureEng = enterEvidence(featureEng, testdata, i);
    [nodePot, edgePot] = mkPotentials(featureEng, weights);
    fprintf([int2str(i) '. ']);
    [nodeBel, MAPlabels, niter, edgeBel, logZ] = infer(infEng, nodePot, edgePot);
    
    if (exist(image_paths_test{i}, 'file'))
        original_image = imread(image_paths_test{i});
        
        [height width channels] = size(original_image);
        image_name = regexprep(image_paths_test{i}, '(.*)/', '');
        
        filenames{i} = image_name;
        
        imwrite(original_image, [output_directory '/' int2str(round(logZ)) '_' image_name]); 
       
        crf_inference_map = imresize(MAPlabels, [height width]);
        crf_inference_map = crf_inference_map / max(crf_inference_map(:));

        %int2str(round(logZ))
        fprintf([output_directory '/' int2str(round(logZ)) '_' image_name '.infer.jpg' '\n']);
        imwrite(crf_inference_map, [output_directory '/' int2str(round(logZ)) '_' image_name '.infer.jpg']);
        
        scores(i) = round(logZ);
    
        already_processed{i} = [output_directory '/' int2str(round(logZ)) '_' image_name '.infer.jpg'];
    end
  
end

save(['saliency_scores_' prefix '.mat'], 'filenames', 'scores');


