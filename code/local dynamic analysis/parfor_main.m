% 定义基本参数结构体 p_base
p_base.tau_N = 0.1; 
p_base.tau_A = 0.002; 
p_base.tau_G = 0.01; 
p_base.tau_r = 0.02;
p_base.alpha = 4.97540390783020; 
p_base.beta = 17.8040082933509; 
p_base.gamma = 0.717362634377862;
p_base.W_E = 1; 
p_base.W_I = 0.6; 
p_base.Ib_E = 0.325; 
p_base.Ib_I = 0.3;
p_base.aE = 310.; 
p_base.bE = 125.; 
p_base.dE = 0.16;
p_base.aI = 615.; 
p_base.bI = 177.; 
p_base.dI = 0.087;

% 计算参数向量
wee_NMDA = 1.46300335335758 * NMDA;
wee_AMPA = 5.42080840436637 * AMPA;
wie = 2.15375172560333 * GABA;
wei_NMDA = 0.335953964409294 * NMDA;
wei_AMPA = 10.8178060053020 * AMPA;
wii = 0.197708592716577 * GABA;

% 确定扫描长度（假设这些向量长度相同）
N = length(wee_NMDA);  % 或者 length(wee_AMPA), 等等

% 预分配存储结果的数组
rE_results = zeros(N, 1);
rI_results = zeros(N, 1);
max_real_eig = zeros(N, 1);
min_real_eig = zeros(N, 1);      % 最小实部
num_complex_pairs = zeros(N, 1); % 复数对的个数
eig_real_parts = cell(N, 1);     % 所有实部
eig_imag_parts = cell(N, 1);     % 所有虚部
equilibrium_type = zeros(N, 1);  % 类型编号
has_converged = false(N, 1);

tic
% 使用 parfor 并行循环
parfor i = 1:N
    % 为每个迭代创建参数结构的副本
    p = p_base;
    
    % 设置当前迭代的参数
    p.omega_N_EE = wee_NMDA(i);
    p.omega_A_EE = wee_AMPA(i);
    p.omega_G_IE = wie(i);
    p.omega_N_EI = wei_NMDA(i);
    p.omega_A_EI = wei_AMPA(i);
    p.omega_G_II = wii(i);
    
    % 初始猜测（可以根据需要调整，比如使用上一个平衡点作为下一个的初值）
    % 但在 parfor 中不能直接依赖前一次结果，所以用固定初值或自适应策略
    r0 = [0.5; 0.5];  % 固定初始猜测
    
    try
        % 求平衡点
        [rE, rI] = find_equilibrium(p, r0);
        
        % 计算雅可比矩阵
        J = compute_jacobian(p, rE, rI);
        
        % 特征值分析
        eigs_val = eig(J);
        real_parts = real(eigs_val);
        imag_parts = imag(eigs_val);

        % 统计信息
        max_real = max(real_parts);
        min_real = min(real_parts);

        % 计算复数对的个数（虚部非零，且只计一次）
        imag_nonzero = abs(imag_parts) > 1e-10;
        num_pairs = sum(imag_nonzero) / 2;  % 应该是整数
        
        % 判断平衡点类型
        type = classify_equilibrium(real_parts, imag_parts);
        
        % 存储结果
        rE_results(i) = rE;
        rI_results(i) = rI;
        max_real_eig(i) = max_real;
        min_real_eig(i) = min_real;
        num_complex_pairs(i) = num_pairs;
        eig_real_parts{i} = real_parts;
        eig_imag_parts{i} = imag_parts;
        equilibrium_type(i) = type;
        has_converged(i) = true;
        
    catch ME
        warning('参数索引 %d 计算失败: %s', i, ME.message);
        rE_results(i) = NaN;
        rI_results(i) = NaN;
        max_real_eig(i) = NaN;
        min_real_eig(i) = NaN;
        num_complex_pairs(i) = NaN;
        eig_real_parts{i} = [];
        eig_imag_parts{i} = [];
        equilibrium_type(i) = 0;  % 0 表示失败
        has_converged(i) = false;
    end
end
toc

% 统计每种类型的数量
converged_idx = has_converged;
valid_types = equilibrium_type(converged_idx);
unique_types = unique(valid_types);
counts = histcounts(valid_types, [unique_types; max(unique_types)+1]);

% 定义分类函数
function type = classify_equilibrium(real_parts, imag_parts)
    % type 编号:
    % 1: 稳定结点
    % 2: 不稳定结点
    % 3: 鞍点
    % 4: 稳定焦点
    % 5: 不稳定焦点
    % 6: 鞍-焦点
    % 7: 其他
    
    tol = 1e-8;
    n = length(real_parts);

    is_real = abs(imag_parts) < tol;
    num_real = sum(is_real);
    num_complex = n - num_real;

    real_eigs = real_parts(is_real);
    complex_real_parts = real_parts(~is_real);

    % 1) 全实特征值
    if num_real == n
        if all(real_eigs < -tol)
            type = 1;   % 稳定结点
        elseif all(real_eigs > tol)
            type = 2;   % 不稳定结点
        elseif any(real_eigs < -tol) && any(real_eigs > tol)
            type = 3;   % 鞍点
        else
            type = 7;   % 有接近0的实特征值
        end
        return;
    end

    % 2) 3个实特征值 + 1对复共轭特征值
    if num_real == 3 && num_complex == 2
        pair_real = mean(complex_real_parts);

        if all(real_eigs < -tol) && pair_real < -tol
            type = 4;   % 稳定焦点
        elseif all(real_eigs > tol) && pair_real > tol
            type = 5;   % 不稳定焦点
        else
            type = 6;   % 鞍-焦点
        end
        return;
    end

    % 3) 1个实特征值 + 2对复共轭特征值
    if num_real == 1 && num_complex == 4
        real_eig = real_eigs(1);

        pair_real_vals = cluster_real_parts(complex_real_parts, 1e-6);

        if length(pair_real_vals) ~= 2
            type = 7;
            return;
        end

        if real_eig < -tol && all(pair_real_vals < -tol)
            type = 4;   % 稳定焦点
        elseif real_eig > tol && all(pair_real_vals > tol)
            type = 5;   % 不稳定焦点
        else
            type = 6;   % 鞍-焦点
        end
        return;
    end

    % 4) 其他情况
    type = 7;
end

function group_vals = cluster_real_parts(vals, tol)
    vals = sort(vals(:));
    group_vals = [];

    for i = 1:length(vals)
        if isempty(group_vals)
            group_vals = vals(i);
        else
            if abs(vals(i) - group_vals(end)) > tol
                group_vals(end+1,1) = vals(i);
            end
        end
    end
end