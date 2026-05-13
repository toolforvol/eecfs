% function [mean_knn_accuracy] = cv_nbknnsvm_EECFS(data, class_label, cv_indices, dataset, alpha, maxK)
% 
% Algorithm = 'EECFS';
% 
% knn_accuracy = zeros(1, 10);
% svm_accuracy = zeros(1, 10);
% nb_accuracy = zeros(1, 10);
% 
% time = zeros(1, 10);
% test = zeros(1, 10);
% 
% mb_file_path = strcat('../../result/Real_world_dataset/', Algorithm, '_MB_', dataset, '.txt');
% fprintf('%s\n\n', mb_file_path);
% fid_mb = fopen(mb_file_path, 'wt');
% 
% 
% for i = 1:10
%     i
%     test_indices = (cv_indices == i);
%     train_indices = ~test_indices;
% 
%     train_data = [data(train_indices, :), class_label(train_indices)];
%     [~, target] = size(train_data);
%     [~, p] = size(train_data);
% 
%     ns = max(train_data);
% 
%     % Retrieve Markov Blanket
%     [mb, ntest, ntime] = EECFS_G2(train_data, target, alpha, ns, p, maxK);
%     time(i) = ntime;
%     test(i) = ntest;
% 
%     lengthMB(i) = length(mb);
%     % Write to file
%     % fprintf(fid_mb, '%s\n', num2str(mb));
%     fprintf(fid_mb, '%s\n', strjoin(string(mb), ' '));
% 
%     if isempty(mb)
%         knn_accuracy(i) = 0;
%         nb_accuracy(i) = 0;
%         svm_accuracy(i) = 0;
%     else
%         % Retrieve class number
%         classes = unique(class_label);
%         num_classes = length(classes);
% 
%         % Dynamically set KNN NumNeighbors parameter
%         knn_model = fitcknn(data(train_indices, mb), class_label(train_indices), 'NumNeighbors', 3);
%         test_class = predict(knn_model, data(test_indices, mb));
%         knn_accuracy(i) = length(find(class_label(test_indices) == test_class)) / length(test_class);
% 
%         if num_classes > 2
%             % For multi-class, use one-vs-rest voting
%             models = cell(num_classes, 1);
%             predicted_labels = cell(num_classes, 1);
% 
%             % Train multiple binary classifiers
%             for m = 1:num_classes
%                 idx = (class_label(train_indices) == classes(m));
%                 data_binary = data(train_indices, mb);
%                 data_binary = data_binary(idx, :);
% 
%                 binary_labels = class_label(train_indices);
%                 binary_labels = binary_labels(idx);
% 
%                 binary_labels(binary_labels == classes(m)) = 1;
%                 binary_labels(binary_labels ~= 1) = -1;
% 
%                 models{m} = fitcsvm(data_binary, binary_labels, 'KernelFunction', 'rbf');
%             end
% 
%             % Predict on test set and perform voting
%             test_data = data(test_indices, mb);
%             test_labels = class_label(test_indices);
% 
%             predicted_label = zeros(sum(test_indices), 1);
%             votes = zeros(sum(test_indices), num_classes);
% 
%             for m = 1:num_classes
%                 test_indices_binary = (test_labels == classes(m));
% 
%                 if any(test_indices_binary)
%                     predicted_labels{m} = zeros(sum(test_indices), 1);
%                     predicted_labels{m}(test_indices_binary) = predict(models{m}, test_data(test_indices_binary, :));
% 
%                     for k = 1:sum(test_indices)
%                         if predicted_labels{m}(k) == 1
%                             votes(k, m) = votes(k, m) + 1;
%                         end
%                     end
%                 end
%             end
% 
%             for k = 1:sum(test_indices)
%                 [~, predicted_class] = max(votes(k, :));
%                 predicted_label(k) = classes(predicted_class);
%             end
% 
%             svm_accuracy(i) = length(find(test_labels == predicted_label)) / length(predicted_label);
% 
%         else
%             % For binary classification, train binary classifier directly
%             model = fitcsvm(data(train_indices, mb), class_label(train_indices), 'KernelFunction', 'rbf');
%             predicted_label = predict(model, data(test_indices, mb));
%             svm_accuracy(i) = length(find(class_label(test_indices) == predicted_label)) / length(predicted_label);
%         end
%     end
% end
% 
% fclose(fid_mb);
% 
% std_knn_accuracy = std(knn_accuracy);
% mean_knn_accuracy = mean(knn_accuracy);
% std_svm_accuracy = std(svm_accuracy);
% mean_svm_accuracy = mean(svm_accuracy);
% std_length = std(lengthMB);
% mean_length = mean(lengthMB);
% std_test = std(test);
% mean_test = mean(test);
% std_time = std(time);
% mean_time = mean(time);
% 
% fprintf('\nKNN\tSVM\tCompactness\tTest\tTime\n\n');
% fprintf('%.2f+%.2f\t%.2f+%.2f\t%.0f+%.0f\t%.0f+%.0f\t%.4f+%.4f\n\n\n', ...
%     mean_knn_accuracy, std_knn_accuracy, mean_svm_accuracy, std_svm_accuracy, ...
%     mean_length, std_length, mean_test, std_test, mean_time, std_time);
% 
% result_str = strcat(strcat('../../result/Real_world_dataset/', Algorithm, '_metrics_', dataset), '.txt');
% 
% fid = fopen(result_str, 'wt');
% fprintf(fid, 'KNN\tSVM\tCompactness\tTest\tTime\n');
% fprintf(fid, '%.2f+%.2f\t%.2f+%.2f\t%.0f+%.0f\t%.0f+%.0f\t%.4f+%.4f', ...
%     mean_knn_accuracy, std_knn_accuracy, mean_svm_accuracy, std_svm_accuracy, ...
%     mean_length, std_length, mean_test, std_test, mean_time, std_time);
% fclose(fid);
% 
% % ++++++++++ 添加CI区间-开始 ++++++++++
% n = length(knn_accuracy);
% t_value = tinv(0.975, n-1);
% ci_knn = t_value * std_knn_accuracy / sqrt(n);
% ci_svm = t_value * std_svm_accuracy / sqrt(n);
% lower_knn = mean_knn_accuracy - ci_knn;
% upper_knn = mean_knn_accuracy + ci_knn;
% lower_svm = mean_svm_accuracy - ci_svm;
% upper_svm = mean_svm_accuracy + ci_svm;
% 
% % Save CI to file
% ci_file = strcat('../../result/Real_world_dataset/', ...
%                  Algorithm, '_CI_', dataset, '.txt');
% fid_ci = fopen(ci_file, 'wt');
% fprintf(fid_ci, ...
% 'Algorithm\tDataset\tClassifier\tMean\tCI\tLower\tUpper\n');
% fprintf(fid_ci, ...
% '%s\t%s\tKNN\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
% Algorithm, dataset, ...
% mean_knn_accuracy, ci_knn, ...
% lower_knn, upper_knn);
% fprintf(fid_ci, ...
% '%s\t%s\tSVM\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
% Algorithm, dataset, ...
% mean_svm_accuracy, ci_svm, ...
% lower_svm, upper_svm);
% fclose(fid_ci);
% % ++++++++++ 添加CI区间-结束 ++++++++++
% 
% 
% end

% =============================[嵌套交叉验证]==============================================

function [mean_knn_accuracy] = cv_nbknnsvm_EECFS( ...
    data, class_label, cv_indices, dataset, alpha, maxK)

Algorithm = 'EECFS';

outerK = 10;

knn_accuracy = zeros(1, outerK);
svm_accuracy = zeros(1, outerK);
nb_accuracy = zeros(1, outerK);

time = zeros(1, outerK);
test = zeros(1, outerK);
lengthMB = zeros(1, outerK);

%% ==========================================================
% Read precomputed MB file
%% ==========================================================
mb_file_path = strcat( ...
    '../../result/Real_world_dataset/', ...
    Algorithm, '_MB_', dataset, '.txt');

fprintf('Load MB from: %s\n\n', mb_file_path);

fid_mb = fopen(mb_file_path, 'rt');

if fid_mb == -1
    error('Cannot open MB file: %s', mb_file_path);
end

mb_list = cell(outerK, 1);

line_idx = 1;
while ~feof(fid_mb)
    line = fgetl(fid_mb);

    if ischar(line) && ~isempty(strtrim(line))
        mb_list{line_idx} = str2num(line); %#ok<ST2NM>

        line_idx = line_idx + 1;
    end
end

fclose(fid_mb);

if length(mb_list) < outerK
    error('MB file contains fewer than %d folds.', outerK);
end

%% ==========================================================
% Outer CV (Nested CV Evaluation)
%% ==========================================================
for i = 1:outerK

    fprintf('Outer Fold %d\n', i);

    test_indices = (cv_indices == i);
    train_indices = ~test_indices;

    %% ------------------------------------------
    % Use precomputed MB for current fold
    %% ------------------------------------------
    mb = mb_list{i};

    % Placeholder (since no rerun EECFS_G2)
    ntime = 0;
    ntest = 0;

    time(i) = ntime;
    test(i) = ntest;

    lengthMB(i) = length(mb);

    %% ------------------------------------------
    % Empty MB
    %% ------------------------------------------
    if isempty(mb)

        knn_accuracy(i) = 0;
        nb_accuracy(i) = 0;
        svm_accuracy(i) = 0;

        continue;
    end

    %% ------------------------------------------
    % Retrieve class number
    %% ------------------------------------------
    classes = unique(class_label);
    num_classes = length(classes);

    %% ======================================================
    % KNN
    %% ======================================================
    knn_model = fitcknn( ...
        data(train_indices, mb), ...
        class_label(train_indices), ...
        'NumNeighbors', 3);

    test_class = predict( ...
        knn_model, ...
        data(test_indices, mb));

    knn_accuracy(i) = ...
        sum(class_label(test_indices) == test_class) ...
        / length(test_class);

    %% ======================================================
    % SVM
    %% ======================================================
    if num_classes > 2

        % One-vs-rest multi-class SVM
        models = cell(num_classes, 1);
        predicted_labels = cell(num_classes, 1);

        for m = 1:num_classes

            binary_labels = class_label(train_indices);

            binary_labels( ...
                binary_labels == classes(m)) = 1;

            binary_labels( ...
                binary_labels ~= 1) = -1;

            models{m} = fitcsvm( ...
                data(train_indices, mb), ...
                binary_labels, ...
                'KernelFunction', 'rbf');
        end

        test_data = data(test_indices, mb);
        test_labels = class_label(test_indices);

        predicted_label = zeros(sum(test_indices), 1);
        votes = zeros(sum(test_indices), num_classes);

        for m = 1:num_classes

            pred = predict( ...
                models{m}, ...
                test_data);

            votes(:, m) = (pred == 1);
        end

        [~, predicted_class] = max(votes, [], 2);

        predicted_label = ...
            classes(predicted_class);

        svm_accuracy(i) = ...
            sum(test_labels == predicted_label) ...
            / length(predicted_label);

    else

        % Binary SVM
        model = fitcsvm( ...
            data(train_indices, mb), ...
            class_label(train_indices), ...
            'KernelFunction', 'rbf');

        predicted_label = predict( ...
            model, ...
            data(test_indices, mb));

        svm_accuracy(i) = ...
            sum(class_label(test_indices) ...
            == predicted_label) ...
            / length(predicted_label);
    end
end

%% ==========================================================
% Statistics
%% ==========================================================
std_knn_accuracy = std(knn_accuracy);
mean_knn_accuracy = mean(knn_accuracy);

std_svm_accuracy = std(svm_accuracy);
mean_svm_accuracy = mean(svm_accuracy);

std_length = std(lengthMB);
mean_length = mean(lengthMB);

std_test = std(test);
mean_test = mean(test);

std_time = std(time);
mean_time = mean(time);

fprintf('\nKNN\tSVM\tCompactness\tTest\tTime\n\n');

fprintf( ...
'%.2f+%.2f\t%.2f+%.2f\t%.0f+%.0f\t%.0f+%.0f\t%.4f+%.4f\n\n\n', ...
mean_knn_accuracy, std_knn_accuracy, ...
mean_svm_accuracy, std_svm_accuracy, ...
mean_length, std_length, ...
mean_test, std_test, ...
mean_time, std_time);

%% ==========================================================
% Save Metrics
%% ==========================================================
result_str = strcat( ...
'../../result/revise/Real_world_dataset/', ...
Algorithm, '_metrics_', dataset, '.txt');

fid = fopen(result_str, 'wt');

fprintf(fid, ...
'KNN\tSVM\tCompactness\tTest\tTime\n');

fprintf(fid, ...
'%.2f+%.2f\t%.2f+%.2f\t%.0f+%.0f\t%.0f+%.0f\t%.4f+%.4f', ...
mean_knn_accuracy, std_knn_accuracy, ...
mean_svm_accuracy, std_svm_accuracy, ...
mean_length, std_length, ...
mean_test, std_test, ...
mean_time, std_time);

fclose(fid);

%% ==========================================================
% 95% CI
%% ==========================================================
n = length(knn_accuracy);

t_value = tinv(0.975, n - 1);

ci_knn = ...
    t_value * std_knn_accuracy / sqrt(n);

ci_svm = ...
    t_value * std_svm_accuracy / sqrt(n);

lower_knn = mean_knn_accuracy - ci_knn;
upper_knn = mean_knn_accuracy + ci_knn;

lower_svm = mean_svm_accuracy - ci_svm;
upper_svm = mean_svm_accuracy + ci_svm;

ci_file = strcat( ...
    '../../result/revise/Real_world_dataset/', ...
    Algorithm, '_CI_', dataset, '.txt');

fid_ci = fopen(ci_file, 'wt');

fprintf(fid_ci, ...
'Algorithm\tDataset\tClassifier\tMean\tCI\tLower\tUpper\n');

fprintf(fid_ci, ...
'%s\t%s\tKNN\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
Algorithm, dataset, ...
mean_knn_accuracy, ci_knn, ...
lower_knn, upper_knn);

fprintf(fid_ci, ...
'%s\t%s\tSVM\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
Algorithm, dataset, ...
mean_svm_accuracy, ci_svm, ...
lower_svm, upper_svm);

fclose(fid_ci);

end
