ans_idx = find(pf > 8 & pf < 13);

if numel(ans_idx) >= 2
    % 分箱并求平均
    loc_y = double(loc(ans_idx, 2));
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
    pf1_clean = pf1(valid);

    % 如果所有点都被剔除，则报错并停止
    if isempty(binIdx_clean)
        error('所有点都落在分箱边界之外，请检查 step 或边界定义。');
    end

    % 4. 计算每个箱的平均值
    loc_mid = (binEdges(1:end-1)+binEdges(2:end))/2;
    % 对没有任何点落入的区间索引，输出用第5个参数nan填充
    pf1Mean = accumarray(binIdx_clean, pf1_clean, [numBins, 1], @mean, NaN);

    valid = find(~isnan(pf1Mean));
    corr_mat = corrcoef(pf1Mean(valid), loc_mid(valid));
    r_pf_y = corr_mat(1, 2);
end
