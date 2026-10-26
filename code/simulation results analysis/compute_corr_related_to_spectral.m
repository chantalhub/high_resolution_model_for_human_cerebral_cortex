% 运行完TotalTimeSeries_SpectrumAnalysis_uptodate_hanningwindow.m，利用以下代码计算相关性
for i = 1:59157
    [~, position] = max(power(i, :));
    pf(i) = frequency(i, position);
end

ratio = AMPA./NMDA;

rng(300);
ans_idx = find(pf > 8 & pf < 13);
indices = randperm(length(ans_idx), 600);
[~, score, ~] = pca(power, 'NumComponents', 10);

if numel(ans_idx) >= 2
    % 分箱并求平均
    loc_y = double(loc(ans_idx, 2));
    pc1 = score(ans_idx, 1);
    t1wt2w1 = t1wt2w;
    t1wt2w1 = t1wt2w1(ans_idx);
    ratio1 = ratio(ans_idx);
    pf1 = pf(ans_idx);
    minY = min(loc_y);
    maxY = max(loc_y);

    % 分箱
    step = 0.02;
    binEdges = minY : step : maxY;

    if binEdges(end) < maxY
        binEdges(end+1) = binEdges(end) + 0.01;
    end
    numBins = length(binEdges) - 1;
    binIdx = discretize(loc_y, binEdges);

    % 3. 移除 binIdx 中为 NaN 的元素（即超出边界的点）
    valid = ~isnan(binIdx);
    
    binIdx_clean = binIdx(valid);
    loc_clean = loc_y(valid);
    pc1_clean = pc1(valid);
    t1wt2w_clean = t1wt2w1(valid);
    ratio1_clean = ratio1(valid);
    pf1_clean = pf1(valid);

    % 如果所有点都被剔除，则报错并停止
    if isempty(binIdx_clean)
        error('所有点都落在分箱边界之外，请检查 step 或边界定义。');
    end

    % 4. 计算每个箱的平均值
    loc_mid = (binEdges(1:end-1)+binEdges(2:end))/2;
    % 对没有任何点落入的区间索引，输出用第5个参数nan填充
    pc1Mean = accumarray(binIdx_clean, pc1_clean, [numBins, 1], @mean, NaN);
    t1wt2wMean = accumarray(binIdx_clean, t1wt2w_clean, [numBins, 1], @mean, NaN);
    ratio1Mean = accumarray(binIdx_clean, ratio1_clean, [numBins, 1], @mean, NaN);
    pf1Mean = accumarray(binIdx_clean, pf1_clean, [numBins, 1], @mean, NaN);

    valid = find(~isnan(pf1Mean));
    corr_mat = corrcoef(-pc1Mean(valid), t1wt2wMean(valid));
    r_t1_score1 = corr_mat(1, 2);
    corr_mat = corrcoef(-pc1Mean(valid), loc_mid(valid));
    r_score1_y = corr_mat(1, 2);
    corr_mat = corrcoef(pf1Mean(valid), ratio1Mean(valid));
    r_pf_ratio = corr_mat(1, 2);
    corr_mat = corrcoef(ratio1Mean(valid), loc_mid(valid));
    r_ratio_y = corr_mat(1, 2);
    corr_mat = corrcoef(pf1Mean(valid), loc_mid(valid));
    r_pf_y = corr_mat(1, 2);
end


