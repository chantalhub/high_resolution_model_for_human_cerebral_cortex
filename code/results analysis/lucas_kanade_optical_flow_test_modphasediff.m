[M, N, P] = size(data);
u_flow = zeros(M, N, P-1);  % 速度x分量
v_flow = zeros(M, N, P-1);  % 速度y分量
flow_direction = zeros(M, N, P-1);  % 传播方向

for t = 1:P-1
    % 将连续两时刻的数据视为"图像帧"
    I1 = data(:,:,t);      % 当前帧
    I2 = data(:,:,t+1);    % 下一帧
    
    % 计算光流
    [u, v, mag, dir] = lucas_kanade_optical_flow(I1, I2, 3);
    
    u_flow(:,:,t) = u;
    v_flow(:,:,t) = v;
    flow_direction(:,:,t) = dir;
end

%{
% 计算整体传播方向（加权平均）
overall_direction = zeros(999, 1);
for t = 1:999
    % 使用速度大小作为权重
    weights = sqrt(u_flow(:,:,t).^2 + v_flow(:,:,t).^2);
    weights = weights / sum(weights(:));  % 归一化
    
    % 加权平均方向
    mean_u = sum(u_flow(:,:,t) .* weights, 'all');
    mean_v = sum(v_flow(:,:,t) .* weights, 'all');
    overall_direction(t) = atan2(mean_v, mean_u);
end
%}

function [u, v, magnitude, direction] = lucas_kanade_optical_flow(I1, I2, window_size)
    % Lucas-Kanade 光流法
    % 输入：I1, I2 - 连续两帧图像，window_size - 窗口大小（如5）
    % 输出：u, v - 速度分量，magnitude - 速度大小，direction - 速度方向（弧度）
    
    % 计算时间梯度
    It = mod(I2 - I1 + pi, 2*pi) - pi;
    
    [Ix, Iy] = gradientMATLAB(I1);
    
    % 初始化速度场
    [height, width] = size(I1);
    u = zeros(height, width);
    v = zeros(height, width);
    direction = NaN(height, width);  % 初始化为NaN
    
    half_window = floor(window_size / 2);
    
    % 对每个像素计算光流
    for i = half_window+1 : height-half_window
        for j = half_window+1 : width-half_window
            % 提取局部窗口
            Ix_window = Ix(i-half_window:i+half_window, j-half_window:j+half_window);
            Iy_window = Iy(i-half_window:i+half_window, j-half_window:j+half_window);
            It_window = It(i-half_window:i+half_window, j-half_window:j+half_window);
            
            % 检测梯度特性
            mean_Ix = mean(abs(Ix_window(:)));
            mean_Iy = mean(abs(Iy_window(:)));
            
            if isnan(mean_Ix) || isnan(mean_Iy)
                u(i,j) = nan;
                v(i,j) = nan;
                
            % 情况1: 两个方向梯度都接近零（无纹理区域）
            elseif mean_Ix < 1e-6 && mean_Iy < 1e-6
                % 梯度为零，无法确定传播方向
                u(i,j) = 0;
                v(i,j) = 0;
            
            % 情况2: y方向梯度接近零，使用一维光流
            elseif mean_Iy < 0.01 * mean_Ix && mean_Ix > 1e-6
                % 一维光流：只考虑x方向
                A = Ix_window(:);
                b = -It_window(:);
                u(i,j) = (A' * A) \ (A' * b);
                v(i,j) = 0;
                
            % 情况3: x方向梯度接近零，使用一维光流  
            elseif mean_Ix < 0.01 * mean_Iy && mean_Iy > 1e-6
                % 一维光流：只考虑y方向
                A = Iy_window(:);
                b = -It_window(:);
                v(i,j) = (A' * A) \ (A' * b);
                u(i,j) = 0;
                
            % 情况4: 正常二维光流
            else
                A = [Ix_window(:), Iy_window(:)];
                b = -It_window(:);
                ATA = A' * A;
                
                if det(ATA) < 1e-10
                    % 如果矩阵奇异，回退到一维方法
                    if mean_Ix > mean_Iy
                        u(i,j) = (Ix_window(:)' * Ix_window(:)) \ (Ix_window(:)' * (-It_window(:)));
                        v(i,j) = 0;
                    else
                        v(i,j) = (Iy_window(:)' * Iy_window(:)) \ (Iy_window(:)' * (-It_window(:)));
                        u(i,j) = 0;
                    end
                else
                    v_vector = ATA \ (A' * b);
                    u(i,j) = v_vector(1);
                    v(i,j) = v_vector(2);
                end
            end
        end
    end
    
    % 计算速度大小和方向
    magnitude = sqrt(u.^2 + v.^2);
    direction = atan2(v, u);
    speed_zero = (abs(u) < 1e-10) & (abs(v) < 1e-10);
    direction(speed_zero) = NaN;
end