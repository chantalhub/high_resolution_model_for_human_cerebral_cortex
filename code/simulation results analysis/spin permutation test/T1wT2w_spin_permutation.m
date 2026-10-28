% 导入/home/liuchang/cb_review/reviewer1/7and13and15/spin_permutation_15中的time_to_peak.mat（record变量）
% 确定左右半球的顶点数 (fs_LR 32k 通常为 32492 每半球)
n_hemi = 32492;
n_total = 2 * n_hemi;  % 64984

% 将数据映射到完整的 64k 空间，缺失顶点填 NaN
Y_left_full  = nan(n_hemi, size(t1wt2w, 2));
Y_right_full = nan(n_hemi, size(t1wt2w, 2));

% idx 是全局索引 (1~64984)，需拆分为左右半球局部索引
% 约定：前 n_hemi 为左半球，后 n_hemi 为右半球
is_left = index <= n_hemi;
is_right = index > n_hemi;
local_idx_left  = index(is_left);
local_idx_right = index(is_right) - n_hemi;  % 转成右半球局部索引 (1~32492)

Y_left_full(local_idx_left, :)  = t1wt2w(is_left, :);
Y_right_full(local_idx_right, :) = t1wt2w(is_right, :);

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
t1wt2w_perm = nan(length(index), n_rep); 
for r = 1:n_rep
    % 合并左右半球为 64984 x T
    Y_full_perm = [Y_rand_cell{1}(:,:,r); Y_rand_cell{2}(:,:,r)];
    % 根据 idx 取出对应顶点
    t1wt2w_perm(:, r) = Y_full_perm(index, :);
end

% 保存置换后的相位时间序列，用于后续计算零分布
%save('spin_permuted_phase_ts.mat', 'phase_ts_perm');

% ===== 背侧通路 dorsal（左半球，181-360）=====
dorsal_left = find(label1 == 182 | label1 == 183 | label1 == 193 | label1 == 196 | label1 == 197 | label1 == 199 | label1 == 200 | label1 == 201 | label1 == 203 | label1 == 209 | label1 == 222 | label1 == 225 | label1 == 226 | label1 == 227 | label1 == 228 | label1 == 229 | label1 == 230 | label1 == 275 | label1 == 297 | label1 == 324 | label1 == 325 | label1 == 326 | label1 == 332 | label1 == 336 | label1 == 337 | label1 == 338 | label1 == 339);

% ===== 背侧通路 dorsal（右半球，1-180）=====
dorsal_right = find(label1 == 2 | label1 == 3 | label1 == 13 | label1 == 16 | label1 == 17 | label1 == 19 | label1 == 20 | label1 == 21 | label1 == 23 | label1 == 29 | label1 == 42 | label1 == 45 | label1 == 46 | label1 == 47 | label1 == 48 | label1 == 49 | label1 == 50 | label1 == 95 | label1 == 117 | label1 == 144 | label1 == 145 | label1 == 146 | label1 == 152 | label1 == 156 | label1 == 157 | label1 == 158 | label1 == 159);

% ===== 腹侧通路 ventral（左半球，181-360）=====
ventral_left = find(label1 == 187 | label1 == 198 | label1 == 202 | label1 == 306 | label1 == 307 | label1 == 312 | label1 == 313 | label1 == 314 | label1 == 315 | label1 == 316 | label1 == 317 | label1 == 318 | label1 == 333 | label1 == 334 | label1 == 335 | label1 == 340 | label1 == 343 | label1 == 357);

% ===== 腹侧通路 ventral（右半球，1-180）=====
ventral_right = find(label1 == 7 | label1 == 18 | label1 == 22 | label1 == 126 | label1 == 127 | label1 == 132 | label1 == 133 | label1 == 134 | label1 == 135 | label1 == 136 | label1 == 137 | label1 == 138 | label1 == 153 | label1 == 154 | label1 == 155 | label1 == 160 | label1 == 163 | label1 == 177);

dorsal = [dorsal_left; dorsal_right];
[r_dorsal, p_dorsal] = corrcoef(t1wt2w(dorsal), record(dorsal)*0.0001);
rho_obs_dorsal = r_dorsal(1,2);
rho_null_dorsal = zeros(n_rep, 1);

for i = 1:n_rep
    t1wt2w_perm1 = t1wt2w_perm(dorsal, i);
    record1 = record(dorsal);
    mask = ~isnan(t1wt2w_perm1);
    t1wt2w_perm2 = t1wt2w_perm1(mask);
    record2 = record1(mask);
    r = corrcoef(t1wt2w_perm2, record2*0.0001);
    rho_null_dorsal(i) = r(1,2);
end
p_spin_dorsal = (sum(rho_null_dorsal<rho_obs_dorsal)+1) / (n_rep+1);

f = figure('Units', 'centimeters', 'Position', [2, 2, 6.4, 4.8]);
histogram(rho_null_dorsal, 30, 'Normalization', 'pdf');
hold on
xline(rho_obs_dorsal, 'r', 'LineWidth', 2);
xlabel('pearson correlation', 'FontName', 'Arial', 'FontSize', 10);
ylabel('probability density', 'FontName', 'Arial', 'FontSize', 10);
%title('dorsal visual processing streams', 'FontName', 'Arial', 'FontSize', 7);

ventral = [ventral_left; ventral_right];
[r_ventral, p_ventral] = corrcoef(t1wt2w(ventral), record(ventral)*0.0001);
rho_obs_ventral = r_ventral(1,2);
rho_null_ventral = zeros(n_rep, 1);

for i = 1:n_rep
    t1wt2w_perm1 = t1wt2w_perm(ventral, i);
    record1 = record(ventral);
    mask = ~isnan(t1wt2w_perm1);
    t1wt2w_perm2 = t1wt2w_perm1(mask);
    record2 = record1(mask);
    r = corrcoef(t1wt2w_perm2, record2*0.0001);
    rho_null_ventral(i) = r(1,2);
end
p_spin_ventral = (sum(rho_null_ventral<rho_obs_ventral)+1) / (n_rep+1);

f = figure('Units', 'centimeters', 'Position', [2, 2, 6.4, 4.8]);
histogram(rho_null_ventral, 30, 'Normalization', 'pdf');
hold on
xline(rho_obs_ventral, 'r', 'LineWidth', 2);
xlabel('pearson correlation', 'FontName', 'Arial', 'FontSize', 10);
ylabel('probability density', 'FontName', 'Arial', 'FontSize', 10);
%title('dorsal visual processing streams', 'FontName', 'Arial', 'FontSize', 7);
