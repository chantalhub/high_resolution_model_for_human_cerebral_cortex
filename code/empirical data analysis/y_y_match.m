T = layout.pos(:,1:2);
X = double(loc(:,1:2));

n = size(X,1);
m = size(T,1);

%% ========= 距离与最近邻 =========
Xn = (X - mean(X,1)) ./ std(X,0,1);
Tn = (T - mean(T,1)) ./ std(T,0,1);

D = pdist2(Xn, Tn, 'squaredeuclidean');
[~, ord] = sort(D, 2, 'ascend');   % 每个点的模板近邻顺序
labels = ord(:,1);                 % 先按纯最近邻分配

counts = accumarray(labels, 1, [m,1]);

% t 为 N×1 向量，labels 为 N×1 向量（取值 1~m）
mean_t1wt2w = accumarray(labels, t1wt2w, [m, 1], @mean, NaN);

non_nan_idx = find(~isnan(mean_t1wt2w));
corr_mat = corrcoef(mean_t1wt2w(non_nan_idx), pc1_scores_mean(non_nan_idx));
r_t1wt2w_score1 = corr_mat(1, 2);
