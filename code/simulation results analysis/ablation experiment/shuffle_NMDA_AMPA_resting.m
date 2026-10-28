% 定义参数和初始条件
I_ext = 0.0;

ind = [1:59157];
nc = length(ind);
left_V1 = evalin('base','left_V1');
right_V1 = evalin('base','right_V1');
input_area = [left_V1; right_V1];

% 所有参数四舍五入保留6位有效数字结果稍逊色
% 只打乱 NMDA/AMPA 比值，不修改 GABA
perm_idx = randperm(nc);
NMDA = NMDA(perm_idx);
AMPA = AMPA(perm_idx);
% GABA 保持不变

wee_NMDA = 1.46300335335758*NMDA; %1.463003*NMDA;
wee_AMPA = 5.42080840436637*AMPA; %5.420808*AMPA;
wie = -2.15375172560333*GABA; %-2.153752*GABA;
wei_NMDA = 0.335953964409294*NMDA; %0.335954*NMDA;
wei_AMPA = 10.8178060053020*AMPA; %10.817806*AMPA;
wii = -0.197708592716577*GABA; %-0.197709*GABA;

sc = evalin('base','sc');

sc = sc ./ sqrt(sum(sc, 2));
sc = sc / mean(sum(sc, 2));

gpuDevice(2);

%euclidean_distance = double(pdist2(loc, loc));
[rows, cols, vals] = find(sc);
vals = gpuArray(vals);
%dist_values = euclidean_distance(sub2ind(size(sc), rows, cols));
dist_values = geodesic_distance(sub2ind(size(sc), rows, cols));

%{
%init_val = evalin('base', 'init_val');
init_val = [0.01 * rand(59157,3), 1 * rand(59157,2)];

S_E_NMDA_init = init_val(ind, 1);
S_E_AMPA_init = init_val(ind, 2);
S_I_init = init_val(ind, 3);
r_E_init = init_val(ind, 4);
r_I_init = init_val(ind, 5);
I_E_init = exc_current(S_I_init, S_E_NMDA_init, S_E_AMPA_init, I_ext, wee_NMDA, wee_AMPA, wie);
I_I_init = inh_current(S_I_init, S_E_NMDA_init, S_E_AMPA_init, wei_NMDA, wei_AMPA, wii);
%I_E_init = init_val(ind, 6);
%I_I_init = init_val(ind, 7);  
%}
init_val = [0.01+0.01*rand(nc,3), 1+1*rand(nc,2), 0.1+0.1*rand(nc,2)];

S_E_NMDA_init = init_val(ind, 1);
S_E_AMPA_init = init_val(ind, 2);
S_I_init = init_val(ind, 3);
r_E_init = init_val(ind, 4);
r_I_init = init_val(ind, 5);
%I_E_init = exc_current(S_I_init, S_E_NMDA_init, S_E_AMPA_init, I_ext, wee_NMDA, wee_AMPA, wie);
%I_I_init = inh_current(S_I_init, S_E_NMDA_init, S_E_AMPA_init, wei_NMDA, wei_AMPA, wii);
I_E_init = init_val(ind, 6);
I_I_init = init_val(ind, 7);

% 设置时间范围（只保留前 4s）
h = 0.0001;
% tspan = [0:0.0001:15];
% n_save = 50000;
tspan = [0:0.0001:4];
n_save = 40000;

% time delay
time_delay = round(dist_values/(20000*h));
%clear euclidean_distance
clear geodesic_distance
time_delay = gpuArray(time_delay);

tic
[nmda, ampa, gaba, ie, ii, re, ri] = Runge_Kutta_four(nc, tspan, h, n_save, S_E_NMDA_init, S_E_AMPA_init, S_I_init, I_E_init, I_I_init, r_E_init, r_I_init, I_ext, wee_NMDA, wee_AMPA, wie, wei_NMDA, wei_AMPA, wii, vals, time_delay, rows, cols, input_area);
toc

function [nmda, ampa, gaba, ie, ii, re, ri] = Runge_Kutta_four(nc, t, h, n_save, s_nmda, s_ampa, si, i_e, i_i, r_e, r_i, I_ext, wee_NMDA, wee_AMPA, wie, wei_NMDA, wei_AMPA, wii, sc_vals, delay, pos_x, pos_y, input_vertices)
len = length(t);

nmda = zeros(nc, n_save);
nmda(:, 1) = s_nmda;
ampa = zeros(nc, n_save);
ampa(:, 1) = s_ampa;
gaba = zeros(nc, n_save);
gaba(:, 1) = si;
ie = zeros(nc, n_save);
ie(:, 1) = gather(i_e);
ii = zeros(nc, n_save);
ii(:, 1) = gather(i_i);
re = zeros(nc, n_save);
re(:, 1) = r_e;
ri = zeros(nc, n_save);
ri(:, 1) = r_i;
ampa_gpu = gpuArray(ampa(:, 1:max(delay)+1));

I_noise = zeros(nc, 1);
sigma_noise = 0.007;
tau_E_AMPA = 0.002;

g_E = 1.8;
g_I = 0.0;

I_ext_before = zeros(nc,1);
I_ext_after = zeros(nc,1); % resting: 不施加外部输入

sigma = 0.01;

for i = 2:len
    I_ext = I_ext_before; % resting: 全程无外部输入
    
    %tic
    pos = 1 + mod((i - 1), n_save);
    
    %{
    for p = 1:2
            for q = 1:2
                if delay(p,q) >= i
                    interaction(p, q) = 0;
                else
                    interaction(p, q) = ampa(p, i-delay(p,q));
                end
            end
        end
    %}
    %disp('step1');
    if i <= max(delay)+1
        %tic
        row = find(delay <= i);
        col = mod(i - 2 - delay(row), max(delay) + 1) + 1; % gpuArray
        %toc
        %tic
        linearIndex = pos_y(row); %pos_x(row) is slow;
        %toc
        %tic
        linearIndex1 = linearIndex + (col - 1) * nc; % gpuArray
        %toc
        %tic
        interaction_vals = ampa_gpu(linearIndex1); % gpuArray
        %toc
        %disp('if-step2');
        %tic
        sc_vals_row = sc_vals(row);
        linearIndex2 = pos_x(row);
        %toc
        %disp('if-step3');
        %tic    
        ampa_input = accumarray(linearIndex2, sc_vals_row .* interaction_vals, [nc, 1], @sum, 0);
        %toc
    else
        %disp('else-step2');
        %tic
        col = mod(i - 2 - delay(row), max(delay) + 1) + 1; % gpuArray
        %toc
        %tic
        linearIndex1 = linearIndex + (col - 1) * nc; % gpuArray
        %toc
        %tic
        interaction_vals = ampa_gpu(linearIndex1); % gpuArray
        %toc
        %disp('else-step3');
        %tic
        ampa_input = accumarray(linearIndex2, sc_vals_row .* interaction_vals, [nc, 1], @sum, 0);
        %toc
    end
    
    pos = 1 + mod((i - 1), n_save);    
    
    I_noise = I_noise * exp(-h / tau_E_AMPA) + sqrt(0.002) / 2 * sigma_noise * sqrt(1 - exp(-(2 * h) / tau_E_AMPA)) * randn(nc, 1);
    
    interaction_e = g_E * ampa_input;
    interaction_i = g_I * ampa_input;
    
    I_E = exc_current(si, s_nmda, s_ampa, I_ext, wee_NMDA, wee_AMPA, wie) + interaction_e + I_noise;
    I_I = inh_current(si, s_nmda, s_ampa, wei_NMDA, wei_AMPA, wii) + interaction_i;
    k1_re = phi_E(r_e, I_E);
    k1_ri = phi_I(r_i, I_I);
    r_E = r_e + 1/2 * h * k1_re;
    r_I = r_i + 1/2 * h * k1_ri;
    k1_nmda = SE_NMDA(nc, s_nmda, r_E);
    k1_ampa = SE_AMPA(nc, s_ampa, r_E);
    k1_gaba = SI(nc, si, r_I);
    NMDA = s_nmda + 1/2 * h * k1_nmda;
    AMPA = s_ampa + 1/2 * h * k1_ampa;
    GABA = si + 1/2 * h * k1_gaba;

    I_E = exc_current(GABA, NMDA, AMPA, I_ext, wee_NMDA, wee_AMPA, wie) + interaction_e + I_noise;
    I_I = inh_current(GABA, NMDA, AMPA, wei_NMDA, wei_AMPA, wii) + interaction_i;
    k2_re = phi_E(r_E, I_E);
    k2_ri = phi_I(r_I, I_I);
    r_E = r_e + 1/2 * h * k2_re;
    r_I = r_i + 1/2 * h * k2_ri;
    k2_nmda = SE_NMDA(nc, NMDA, r_E);
    k2_ampa = SE_AMPA(nc, AMPA, r_E);
    k2_gaba = SI(nc, GABA, r_I);
    NMDA = s_nmda + 1/2 * h * k2_nmda;
    AMPA = s_ampa + 1/2 * h * k2_ampa;
    GABA = si + 1/2 * h * k2_gaba;
    
    I_E = exc_current(GABA, NMDA, AMPA, I_ext, wee_NMDA, wee_AMPA, wie) + interaction_e + I_noise;
    I_I = inh_current(GABA, NMDA, AMPA, wei_NMDA, wei_AMPA, wii) + interaction_i;
    k3_re = phi_E(r_E, I_E);
    k3_ri = phi_I(r_I, I_I);
    r_E = r_e + h * k3_re;
    r_I = r_i + h * k3_ri;
    k3_nmda = SE_NMDA(nc, NMDA, r_E);
    k3_ampa = SE_AMPA(nc, AMPA, r_E);
    k3_gaba = SI(nc, GABA, r_I);
    NMDA = s_nmda + h * k3_nmda;
    AMPA = s_ampa + h * k3_ampa;
    GABA = si + h * k3_gaba;
    
    I_E = exc_current(GABA, NMDA, AMPA, I_ext, wee_NMDA, wee_AMPA, wie) + interaction_e + I_noise;
    I_I = inh_current(GABA, NMDA, AMPA, wei_NMDA, wei_AMPA, wii) + interaction_i;
    ie(:, pos) = gather(I_E);
    ii(:, pos) = gather(I_I);
    k4_re = phi_E(r_E, I_E);
    k4_ri = phi_I(r_I, I_I);
    r_E = r_e + h * k4_re;
    r_I = r_i + h * k4_ri;
    re(:, pos) = gather(r_e + 1/6 * h * (k1_re + 2 * k2_re + 2 * k3_re + k4_re));
    ri(:, pos) = gather(r_i + 1/6 * h * (k1_ri + 2 * k2_ri + 2 * k3_ri + k4_ri));
    k4_nmda = SE_NMDA(nc, NMDA, r_E);
    k4_ampa = SE_AMPA(nc, AMPA, r_E);
    k4_gaba = SI(nc, GABA, r_I);
    nmda(:, pos) = gather(s_nmda + 1/6 * h * (k1_nmda + 2 * k2_nmda + 2 * k3_nmda + k4_nmda));
    ampa(:, pos) = gather(s_ampa + 1/6 * h * (k1_ampa + 2 * k2_ampa + 2 * k3_ampa + k4_ampa));
    gaba(:, pos) = gather(si + 1/6 * h * (k1_gaba + 2 * k2_gaba + 2 * k3_gaba + k4_gaba)); 
        
    s_nmda = nmda(:, pos);
    s_ampa = ampa(:, pos);
    si = gaba(:, pos);
    i_e = ie(:, pos);
    i_i = ii(:, pos);
    r_e = re(:, pos);
    r_i = ri(:, pos);

    %disp('step4');
    %tic
    ampa_gpu(:, mod(i-2-max(delay),max(delay)+1)+1) = s_ampa;
    %toc
    
    if max(isnan(r_e))==1
        'nan error'
    end

    if max(isinf(r_e))==1
        'error'
    end
    %toc % Elapsed time is 0.047896 seconds.
    
    %{
    if ~mod(i, 5000)
        r_E_min = min(min(re(:, pos)));
        r_E_max = max(max(re(:, pos)));
        mean_val = mean(re(:, pos));
        disp(['i = ', num2str(i), ', min(r_E) = ', num2str(r_E_min), ', max(r_E) = ', num2str(r_E_max), ', mean(r_E) = ', num2str(mean_val)]);
    end
    %}
    if ~mod(i, 1000)
        r_E_min = min(mean(re(:, pos-999:pos), 2));
        r_E_mean = mean(mean(re(:, pos-999:pos), 2));
        r_E_max = max(mean(re(:, pos-999:pos), 2));
        r_I_min = min(mean(ri(:, pos-999:pos), 2));
        r_I_mean = mean(mean(ri(:, pos-999:pos), 2));
        r_I_max = max(mean(ri(:, pos-999:pos), 2));
        disp(['i = ', num2str(i), ', min(r_E) = ', num2str(r_E_min), ', mean(r_E) = ', num2str(r_E_mean), ', max(r_E) = ', num2str(r_E_max)]);
        disp(['i = ', num2str(i), ', min(r_I) = ', num2str(r_I_min), ', mean(r_I) = ', num2str(r_I_mean), ', max(r_I) = ', num2str(r_I_max)]);
        %[mean(interaction),median(interaction),max(interaction)]
        s_ampa_max_mean =  mean(max(ampa(:, pos-999:pos), [], 2));
        s_ampa_max_median =  median(max(ampa(:, pos-999:pos), [], 2));
        s_ampa_max_max =  max(max(ampa(:, pos-999:pos), [], 2));
        disp(['i = ', num2str(i), ', mean(s_ampa_max_mean) = ', num2str(s_ampa_max_mean), ', median(s_ampa_max_median) = ', num2str(s_ampa_max_median), ', max(s_ampa_max_median) = ', num2str(s_ampa_max_max)]);
        disp(['i = ', num2str(i), ', mean(interaction_e) = ', num2str(mean(interaction_e)), ', meadian(interaction_e) = ', num2str(median(interaction_e)), ', max(interaction_e) = ', num2str(max(interaction_e))]);
    end
    
    %{
    if ~mod(i, n_save)
        %{
        file_name = ['gating_NMDA_', num2str(i), '.mat'];
        save(file_name, 'nmda');
        file_name = ['gating_AMPA_', num2str(i), '.mat'];
        save(file_name, 'ampa');
        file_name = ['gating_GABA_', num2str(i), '.mat'];
        save(file_name, 'gaba');
        file_name = ['rate_I_', num2str(i), '.mat'];
        save(file_name, 'ri');
        %}
        file_name = ['./rateE/rate_E_', num2str(i), '.mat'];
        save(file_name, 're');
    end
    %}
    
end
end

function out = exc_current(S_I, S_E_NMDA, S_E_AMPA, I_ext, wee_NMDA, wee_AMPA, wie)
    W_E = 1.0;
    I_b = 0.325;
    I0_E = W_E * I_b;
    
    out = I0_E + I_ext + (wee_NMDA .* S_E_NMDA) + (wee_AMPA .* S_E_AMPA) + (wie .* S_I);
end

function out = phi_E(r_E, I_E)
    tau_rE = 0.02;
    
    a_E = 310.;
    b_E = 125.;
    d_E = 0.16;   

    out = -(r_E / tau_rE) + (a_E * I_E - b_E) ./ (1 - exp(-d_E * (a_E * I_E - b_E))) / tau_rE;
    
    %{
    C_E = 310.;
    IL_E = 125.;
    d_E = 0.16;
    rmax_E = 80;

    out = -(r_E / tau_rE) + (C_E * I_E - IL_E) ./ (1 - exp(-d_E * (C_E * I_E - IL_E)) + (C_E * I_E - IL_E) / rmax_E) / tau_rE;
    %}
end

function out = inh_current(S_I, S_E_NMDA, S_E_AMPA, wei_NMDA, wei_AMPA, wii)
    W_I = 0.6;
    I_b = 0.3;
    I0_I = W_I * I_b;

    out = I0_I + (wei_NMDA .* S_E_NMDA) + (wei_AMPA .* S_E_AMPA) + (wii .* S_I);
end

function out = phi_I(r_I, I_I)
    tau_rI = 0.02;
    
    a_I = 615.;
    b_I = 177.;
    d_I = 0.087;

    out = -(r_I / tau_rI) + (a_I * I_I - b_I) ./ (1 - exp(-d_I * (a_I * I_I - b_I))) / tau_rI;
    
    %{
    C_I = 615.;
    IL_I = 177.;
    d_I = 0.087;
    rmax_I = 120;

    out = -(r_I / tau_rI) + (C_I * I_I - IL_I) ./ (1 - exp(-d_I * (C_I * I_I - IL_I)) + (C_I * I_I - IL_I) / rmax_I) / tau_rI;
    %}
end

function out = SE_NMDA(nc, S_E_NMDA, r_E)
    tau_E_NMDA = 0.1;
    %sigma = 0.01;
    x = 4.97540390783020; %4.975404;
    
    out = -(S_E_NMDA / tau_E_NMDA) + (x * r_E) .* (1.0 - S_E_NMDA); %+ sigma * randn(nc, 1);
end

function out = SE_AMPA(nc, S_E_AMPA, r_E)
    tau_E_AMPA = 0.002;
    %sigma = 0.01;
    y = 17.8040082933509; %17.804008;
    
    out = -(S_E_AMPA / tau_E_AMPA) + (y * r_E) .* (1.0 - S_E_AMPA); %+ sigma * randn(nc, 1);
end

function out = SI(nc, S_I, r_I)
    tau_I = 0.01;
    %sigma = 0.01;
    z = 0.717362634377862; %0.717363;

    %out =  -(S_I / tau_I) + (z * r_I) .* (1.0 - S_I) + sigma * randn(nc, 1);
    out =  -(S_I / tau_I) + z * r_I; %+ sigma * randn(nc, 1);
end
