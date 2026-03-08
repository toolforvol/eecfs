clc;
clear;

% === Root directory containing all datasets
root_folder = './';

% Get all items in the root directory
folders = dir(root_folder);

% Keep only directories
folders = folders([folders.isdir]);

% Remove '.' and '..' system directories
folders = folders(~ismember({folders.name}, {'.','..'}));

for i = 1:length(folders)
    
    % Get the current subfolder name and path
    subfolder_name = folders(i).name;
    subfolder_path = fullfile(root_folder, subfolder_name);
    
    fprintf("Processing folder: %s\n", subfolder_name);
    
    % === Path to data_labels.mat
    data_file = fullfile(subfolder_path, 'data_labels.mat');
    
    % Skip this folder if the file does not exist
    if ~exist(data_file, 'file')
        fprintf("  data_labels.mat not found, skipping.\n");
        continue
    end
    
    % === Load the dataset
    data_labels = load(data_file);
    
    % Adjust data and labels (shift by +1 if required by the algorithm)
    data = data_labels.data + 1;
    labels = data_labels.labels + 1;
    
    % === Save the processed data and labels back to the same folder
    save(fullfile(subfolder_path, 'data_labels.mat'), 'data', 'labels');
    
    % === Create 10-fold cross-validation indices
    indices = crossvalind('Kfold', labels, 10);
    
    % === Save the cross-validation indices
    save(fullfile(subfolder_path, 'cv10_indices.mat'), 'indices');
    
end

% Print completion message
fprintf("All datasets processed.\n");