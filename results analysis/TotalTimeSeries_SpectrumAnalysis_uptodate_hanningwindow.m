firing_rate = re(:, 30000:40000);
frequency = zeros(59157, 2001);
power = zeros(59157, 2001);
mean_val = mean(firing_rate,2);
rate = (firing_rate - mean_val) ;
for i = 1:59157
    [f, p] = test_fft(rate(i, :));
    frequency(i, :) = f;
    power(i, :) = p;
end

function [f,power] = test_fft(s)
    s = s';
    dt = 0.0001; %unit:s
    fs = 1/dt; %unit:Hz
    t = 1:length(s);
    frei = 40;
    L = length(t);
    
    % 加入汉宁窗并进行能量补偿
    hann_window = hanning(L);  % 生成汉宁窗
    window_power = sum(hann_window.^2);  % 计算窗函数的能量
    s_windowed = s .* hann_window;  % 应用汉宁窗
    
    % 计算FFT
    NFFT = L * frei;
    Y = fft(s_windowed, NFFT) / (L * frei);
    
    % 能量补偿
    Y = Y * sqrt(L / window_power);
    
    f0 = 1 / (dt * L * frei); % 基频
    f = [0:ceil((NFFT-1)/2)] * f0; % 频率序列
    
    % 计算单边频谱
    P1 = abs(Y(1:ceil((NFFT-1)/2)+1)); % 只取正频率部分
    P1(2:end) = 2 * P1(2:end); % 乘以2（除直流分量外）
    
    % 确保输出大小正确
    if length(P1) > 2001
        power = P1(2:2002); % 取前2001个点（排除直流分量）
        f = f(2:2002);
    end
    
    % 确保输出是行向量
    power = power(:)';
    f = f(:)';
end
