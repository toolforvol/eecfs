% % 示例运行时间矩阵：每行对应一个数据集，每列是一个算法（EDMB, EEDMB）
% % 替换成你真实的实验数据
% % 示例：19个数据集上两个算法的运行时间（秒）
% runtimes = [
%     0.07, 0.06;
%     0.05, 0.03;
%     0.22, 0.06;
%     0.04, 0.03;
%     0.11, 0.05;
%     0.1, 0.08;
%     8.33, 2.35;
%     254.48, 61.49;
%     0.18, 0.2;
%     974.73, 180.55;
%     0.82, 0.59;
%     14.4, 5.95;
%     85.45, 6.87;
%     219.84, 196.77;
%     12.49, 7.54;
%     8640, 8315.15;
%     7105.88, 8121.89;
%     123.78, 63.21;
%     6345.06, 8640
% ];
% 
% % Friedman 检验
% [p, tbl, stats] = friedman(runtimes, 1, 'off');  % 'off' 关闭箱线图绘制
% fprintf('Friedman 检验 p 值: %.4f\n', p);
% 
% % 判断是否拒绝原假设（p < 0.05）
% alpha = 0.05;
% if p < alpha
%     fprintf('→ 拒绝原假设，算法之间存在显著差异，执行 Nemenyi 检验\n');
% 
%     % Nemenyi 事后检验
%     % 获取排名（数值越小表示运行时间越快）
%     ranks = tiedrank(runtimes')';  % 每行是一个数据集，对列算法排序
%     avg_ranks = mean(ranks);
%     k = size(runtimes, 2);  % 算法数
%     N = size(runtimes, 1);  % 数据集数
% 
%     % 计算临界差值 (CD, critical difference) for Nemenyi test
%     q_alpha = 2.326;  % 对于 alpha = 0.05，算法数k=2时，q_alpha=2.326 (查Nemenyi表)
%     CD = q_alpha * sqrt(k*(k+1)/(6*N));
%     fprintf('平均秩: EDMB = %.3f, EEDMB = %.3f\n', avg_ranks(1), avg_ranks(2));
%     fprintf('Nemenyi 临界差值 (CD): %.4f\n', CD);
% 
%     % 比较两个算法之间的秩差
%     diff = abs(avg_ranks(1) - avg_ranks(2));
%     if diff > CD
%         fprintf('→ 平均秩差 %.4f > CD，EDMB 与 EEDMB 存在显著差异。\n', diff);
%     else
%         fprintf('→ 平均秩差 %.4f ≤ CD，EDMB 与 EEDMB 差异不显著。\n', diff);
%     end
% else
%     fprintf('→ 无法拒绝原假设，EDMB 与 EEDMB 在运行时间上无显著差异。\n');
% end

% 示例运行时间矩阵：每行对应一个数据集，每列是一个算法（EDMB, EEDMB）
% 示例：17个数据集上两个算法的运行时间（秒）
runtimes = [
    0.07, 0.06;
    0.05, 0.03;
    0.22, 0.06;
    0.04, 0.03;
    0.11, 0.05;
    8.33, 2.35;
    254.48, 61.49;
    0.18, 0.2;
    974.73, 180.55;
    %604800,39311;
    0.82, 0.59;
    14.4, 5.95;
    85.45, 6.87;
    219.84, 196.77;
    8640, 8315.15;
    7105.88, 8121.89;
    123.78, 63.21
];

% Friedman 检验
[p, tbl, stats] = friedman(runtimes, 1, 'off');  % 'off' 关闭箱线图绘制
fprintf('Friedman 检验 p 值: %.4f\n', p);

alpha = 0.05;
if p < alpha
    fprintf('→ 拒绝原假设，算法之间存在显著差异，执行 Nemenyi 检验\n');

    % === 修改部分：运行时间越小越好 → 对其取负值，再排名，越快排名越高 ===
    inv_runtimes = -runtimes;  % 取负值：小值变大
    ranks = tiedrank(inv_runtimes')';  % 每行是一个数据集，对列算法排序
    avg_ranks = mean(ranks);

    k = size(runtimes, 2);  % 算法数
    N = size(runtimes, 1);  % 数据集数
    
    q_alpha = 2.326;  % 对于 alpha = 0.05，k=2
    CD = q_alpha * sqrt(k*(k+1)/(6*N));
    
    fprintf('平均秩（rank 值越大越好）: EDMB = %.3f, EEDMB = %.3f\n', avg_ranks(1), avg_ranks(2));
    fprintf('Nemenyi 临界差值 (CD): %.4f\n', CD);
    
    diff = abs(avg_ranks(1) - avg_ranks(2));
    if diff > CD
        fprintf('→ 平均秩差 %.4f > CD，EDMB 与 EEDMB 存在显著差异。\n', diff);
    else
        fprintf('→ 平均秩差 %.4f ≤ CD，EDMB 与 EEDMB 差异不显著。\n', diff);
    end
else
    fprintf('→ 无法拒绝原假设，EDMB 与 EEDMB 在运行时间上无显著差异。\n');
end

