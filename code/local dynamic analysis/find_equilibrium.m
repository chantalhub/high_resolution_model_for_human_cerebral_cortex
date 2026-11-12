function [rE_star, rI_star] = find_equilibrium(p, r0)
    % p 是包含所有参数的结构体
    % r0 = [rE0; rI0] 初始猜测
    
    % 定义平衡点方程 F(r) = 0
    F = @(r) equilibrium_equations(r, p);
    
    % 用 fsolve 求解
    opts = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'trust-region-dogleg');
    [r_star, ~, exitflag] = fsolve(F, r0, opts);
    
    if exitflag <= 0
        warning('fsolve 未收敛到平衡点');
    end
    
    rE_star = r_star(1);
    rI_star = r_star(2);
end

function F = equilibrium_equations(r, p)
    rE = r(1);
    rI = r(2);
    
    % 稳态 S
    S_N = (p.alpha * p.tau_N * rE) / (1 + p.alpha * p.tau_N * rE);
    S_A = (p.beta * p.tau_A * rE) / (1 + p.beta * p.tau_A * rE);
    S_G = p.gamma * p.tau_G * rI;
    
    % 电流
    I_E = p.omega_N_EE * S_N + p.omega_A_EE * S_A - p.omega_G_IE * S_G + p.W_E * p.Ib_E;
    I_I = p.omega_N_EI * S_N + p.omega_A_EI * S_A - p.omega_G_II * S_G + p.W_I * p.Ib_I;
    
    % 增益函数
    phi_E = gain_func(I_E, p.aE, p.bE, p.dE);
    phi_I = gain_func(I_I, p.aI, p.bI, p.dI);
    
    % 平衡点方程
    F = [rE - phi_E;
         rI - phi_I];
end

function y = gain_func(I, a, b, d)
    x = a * I - b;
    % 防止分母为零
    if abs(x) < 1e-12
        y = 1 / d;  % 极限值，由洛必达法则得到
    else
        y = x / (1 - exp(-d * x));
    end
end