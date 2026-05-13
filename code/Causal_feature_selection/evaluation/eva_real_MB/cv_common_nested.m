
function [mean_knn_accuracy] = cv_common_nested( ...
    data, class_label, cv_indices, dataset, alpha, maxK)

Algorithm = 'EECFS';
outerK = 10;

knn_accuracy = zeros(1, outerK);
svm_accuracy = zeros(1, outerK);
nb_accuracy = zeros(1, outerK);

time = zeros(1, outerK);
test = zeros(1, outerK);
lengthMB = zeros(1, outerK);

% Default return value
mean_knn_accuracy = NaN;

%% ==========================================================
% Read precomputed MB file
%% ==========================================================
mb_file_path = strcat( ...
    '../../result/Real_world_dataset/', ...
    Algorithm, '_MB_', dataset, '.txt');

fprintf('Load MB from: %s\n\n', mb_file_path);

% ---------- check file existence ----------
if ~exist(mb_file_path, 'file')

    warning( ...
        '[%s] MB file does not exist: %s', ...
        dataset, mb_file_path);

    return;
end

fid_mb = fopen(mb_file_path, 'rt');

if fid_mb == -1

    warning( ...
        '[%s] Cannot open MB file: %s', ...
        dataset, mb_file_path);

    return;
end

mb_list = cell(outerK, 1);
line_idx = 1;

while ~feof(fid_mb)

    line = fgetl(fid_mb);

    if ischar(line) && ~isempty(strtrim(line))

        mb_list{line_idx} = str2num(line); %#ok<ST2NM>
        line_idx = line_idx + 1;

        % Avoid reading >10 folds
        if line_idx > outerK + 1
            break;
        end
    end
end

fclose(fid_mb);

% ---------- check fold number ----------
valid_fold_num = line_idx - 1;

if valid_fold_num < outerK

    warning( ...
        '[%s] MB file only contains %d folds (< %d): %s', ...
        dataset, ...
        valid_fold_num, ...
        outerK, ...
        mb_file_path);

    return;
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

    % Since feature selection already finished
    ntime = 0;
    ntest = 0;

    time(i) = ntime;
    test(i) = ntest;

    lengthMB(i) = length(mb);

    %% ------------------------------------------
    % Empty MB
    %% ------------------------------------------
    if isempty(mb)

        warning( ...
            '[%s] Fold %d has empty MB.', ...
            dataset, i);

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
    try
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

    catch ME

        warning( ...
            '[%s] Fold %d KNN failed: %s', ...
            dataset, i, ME.message);

        knn_accuracy(i) = NaN;
    end

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
std_knn_accuracy = std(knn_accuracy, 'omitnan');
mean_knn_accuracy = mean(knn_accuracy, 'omitnan');

std_svm_accuracy = std(svm_accuracy, 'omitnan');
mean_svm_accuracy = mean(svm_accuracy, 'omitnan');

std_length = std(lengthMB, 'omitnan');
mean_length = mean(lengthMB, 'omitnan');

std_test = std(test, 'omitnan');
mean_test = mean(test, 'omitnan');

std_time = std(time, 'omitnan');
mean_time = mean(time, 'omitnan');

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

if fid ~= -1

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

else

    warning( ...
        '[%s] Cannot save metrics file.', ...
        dataset);
end

%% ==========================================================
% 95% Confidence Interval
%% ==========================================================
valid_knn = knn_accuracy(~isnan(knn_accuracy));
valid_svm = svm_accuracy(~isnan(svm_accuracy));

n_knn = length(valid_knn);
n_svm = length(valid_svm);

if n_knn > 1

    t_knn = tinv(0.975, n_knn - 1);

    ci_knn = ...
        t_knn * std(valid_knn) / sqrt(n_knn);

    lower_knn = mean_knn_accuracy - ci_knn;
    upper_knn = mean_knn_accuracy + ci_knn;

else

    ci_knn = NaN;
    lower_knn = NaN;
    upper_knn = NaN;
end

if n_svm > 1

    t_svm = tinv(0.975, n_svm - 1);

    ci_svm = ...
        t_svm * std(valid_svm) / sqrt(n_svm);

    lower_svm = mean_svm_accuracy - ci_svm;
    upper_svm = mean_svm_accuracy + ci_svm;

else

    ci_svm = NaN;
    lower_svm = NaN;
    upper_svm = NaN;
end

%% ==========================================================
% Save CI
%% ==========================================================
ci_file = strcat( ...
    '../../result/revise/Real_world_dataset/', ...
    Algorithm, '_CI_', dataset, '.txt');

fid_ci = fopen(ci_file, 'wt');

if fid_ci ~= -1

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

else

    warning( ...
        '[%s] Cannot save CI file.', ...
        dataset);
end

end