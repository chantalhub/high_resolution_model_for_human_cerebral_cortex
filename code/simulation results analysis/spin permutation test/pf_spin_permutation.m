% 确定左右半球的顶点数 (fs_LR 32k 通常为 32492 每半球)
n_hemi = 32492;
n_total = 2 * n_hemi;  % 64984

% 将数据映射到完整的 64k 空间，缺失顶点填 NaN
pf = pf';
Y_left_full  = nan(n_hemi, size(pf, 2));
Y_right_full = nan(n_hemi, size(pf, 2));

% idx 是全局索引 (1~64984)，需拆分为左右半球局部索引
% 约定：前 n_hemi 为左半球，后 n_hemi 为右半球
is_left = index <= n_hemi;
is_right = index > n_hemi;
local_idx_left  = index(is_left);
local_idx_right = index(is_right) - n_hemi;  % 转成右半球局部索引 (1~32492)

Y_left_full(local_idx_left, :)  = pf(is_left, :);
Y_right_full(local_idx_right, :) = pf(is_right, :);

% 读取球面文件
L_sphere = '/home/liuchang/cb_review/reviewer1/7and13and15/Glasser_et_al_2016_HCP_MMP1.0_v6_RVVG/Glasser_et_al_2016_HCP_MMP1.0_v6_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedValidation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedValidation210.L.sphere.32k_fs_LR.surf.gii';
R_sphere = '/home/liuchang/cb_review/reviewer1/7and13and15/Glasser_et_al_2016_HCP_MMP1.0_v6_RVVG/Glasser_et_al_2016_HCP_MMP1.0_v6_RVVG/Q1-Q6_RelatedValidation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedValidation210.R.sphere.32k_fs_LR.surf.gii';

% 设定旋转次数
n_rep = 1000;

% 调用 spin_permutations（输入为 cell，左右分开）
Y_rand_cell = spin_permutations({Y_left_full, Y_right_full}, ...
                                {L_sphere, R_sphere}, n_rep, ...
                                'surface_algorithm', 'freesurfer');
% 输出：Y_rand_cell{1}: 左半球 32492 x T x 1000
%       Y_rand_cell{2}: 右半球 32492 x T x 1000

% 提取原始59157顶点的置换后数据，并合并半球
pf_perm = nan(length(index), n_rep); 
for r = 1:n_rep
    % 合并左右半球为 64984 x T
    Y_full_perm = [Y_rand_cell{1}(:,:,r); Y_rand_cell{2}(:,:,r)];
    % 根据 idx 取出对应顶点
    pf_perm(:, r) = Y_full_perm(index, :);
end

% 保存置换后的相位时间序列，用于后续计算零分布
%save('spin_permuted_phase_ts.mat', 'phase_ts_perm');

ans = find(pf>8 & pf<13);
r1 = corrcoef(loc(ans, 2), pf(ans));
ratio = AMPA./NMDA;
r2 = corrcoef(ratio(ans), pf(ans));
rho_obs_all1 = r1(1,2);
rho_null_all1 = zeros(n_rep, 1);
for i = 1:n_rep
    mask = ~isnan(pf_perm(:, i));
    pf_perm1 = pf_perm(mask, i);
    loc1 = loc(mask, 2);
    ans1 = find(pf_perm1>8 & pf_perm1<13);
    r = corrcoef(pf_perm1(ans1), loc1(ans1));
    rho_null_all1(i) = r(1,2);
end
p_spin_all1 = (sum(rho_null_all1<rho_obs_all1)+1) / (n_rep+1);

f = figure('Units', 'centimeters', 'Position', [2, 2, 6.4, 4.8]);
histogram(rho_null_all1, 30, 'Normalization', 'pdf');
hold on
xline(rho_obs_all1, 'r', 'LineWidth', 2);
xlabel('spearman correlation');
ylabel('probability density');
set(gca, 'FontName', 'Arial', 'FontSize', 10);

rho_obs_all2 = r2(1,2);
rho_null_all2 = zeros(n_rep, 1);
for i = 1:n_rep
    mask = ~isnan(pf_perm(:, i));
    pf_perm1 = pf_perm(mask, i);
    ratio1 = ratio(mask);
    ans1 = find(pf_perm1>8 & pf_perm1<13);
    r = corrcoef(pf_perm1(ans1), ratio1(ans1));
    rho_null_all2(i) = r(1,2);
end
p_spin_all2 = (sum(rho_null_all2<rho_obs_all2)+1) / (n_rep+1);

f = figure('Units', 'centimeters', 'Position', [2, 2, 6.4, 4.8]);
histogram(rho_null_all2, 30, 'Normalization', 'pdf');
hold on
xline(rho_obs_all2, 'r', 'LineWidth', 2);
xlabel('spearman correlation');
ylabel('probability density');
set(gca, 'FontName', 'Arial', 'FontSize', 10);