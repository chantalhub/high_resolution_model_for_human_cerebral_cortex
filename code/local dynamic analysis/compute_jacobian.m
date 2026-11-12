function J = compute_jacobian(p, rE, rI)
    % 先计算平衡点处的稳态 S
    S_N = (p.alpha * p.tau_N * rE) / (1 + p.alpha * p.tau_N * rE);
    S_A = (p.beta * p.tau_A * rE) / (1 + p.beta * p.tau_A * rE);
    S_G = p.gamma * p.tau_G * rI;
    
    % 计算电流
    I_E = p.omega_N_EE * S_N + p.omega_A_EE * S_A - p.omega_G_IE * S_G + p.W_E * p.Ib_E;
    I_I = p.omega_N_EI * S_N + p.omega_A_EI * S_A - p.omega_G_II * S_G + p.W_I * p.Ib_I;
    
    % 增益函数及其导数
    [phi_E, dphi_E] = gain_func_deriv(I_E, p.aE, p.bE, p.dE);
    [phi_I, dphi_I] = gain_func_deriv(I_I, p.aI, p.bI, p.dI);
    
    % 雅可比矩阵 J 5x5，变量顺序 [S_N, S_A, S_G, rE, rI]
    J = zeros(5);
    
    % 第1行：dS_N/dt
    J(1,1) = -1/p.tau_N - p.alpha * rE;
    J(1,4) = p.alpha * (1 - S_N);
    
    % 第2行：dS_A/dt
    J(2,2) = -1/p.tau_A - p.beta * rE;
    J(2,4) = p.beta * (1 - S_A);
    
    % 第3行：dS_G/dt
    J(3,3) = -1/p.tau_G;
    J(3,5) = p.gamma;
    
    % 第4行：drE/dt (乘以 tau_r 放在右边)
    J(4,1) = dphi_E * p.omega_N_EE / p.tau_r;
    J(4,2) = dphi_E * p.omega_A_EE / p.tau_r;
    J(4,3) = -dphi_E * p.omega_G_IE / p.tau_r;
    J(4,4) = -1/p.tau_r;
    
    % 第5行：drI/dt
    J(5,1) = dphi_I * p.omega_N_EI / p.tau_r;
    J(5,2) = dphi_I * p.omega_A_EI / p.tau_r;
    J(5,3) = -dphi_I * p.omega_G_II / p.tau_r;
    J(5,5) = -1/p.tau_r;
end

function [phi, dphi] = gain_func_deriv(I, a, b, d)
    x = a * I - b;
    if abs(x) < 1e-12
        % x -> 0 时，使用洛必达法则推导的精确极限
        phi = 1/d;
        dphi_dx = 1/2; % 导数的极限：dφ/dx -> 1/2
        dphi = a * dphi_dx; % 链式法则：dφ/dI = a * dφ/dx
    else
        expdx = exp(-d * x);
        phi = x / (1 - expdx);
        % 解析导数
        dphi_dx = (1 - expdx - d * x * expdx) / (1 - expdx)^2;
        dphi = a * dphi_dx;
    end
end