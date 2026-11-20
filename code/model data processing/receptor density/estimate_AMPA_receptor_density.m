load('data\t1wt2w.mat', 't1wt2w');

load('data\loc_inflated.mat', 'loc');
loc = double(loc);

load('data\NMDA.mat', 'NMDA');
load('data\GABA.mat', 'GABA');
NMDA = double(NMDA);
GABA = double(GABA);

NMDA = sqrt(NMDA);
GABA = sqrt(GABA);

y1 = -102.2; % vertex 21752
t1wt2w1 =  1.55;
y2 = 66.5; % vertex 55270
t1wt2w2 =  1.18;
k_slope = (t1wt2w2 - t1wt2w1) / (y2 - y1);
b_intercept = y1 - k_slope * t1wt2w1;
y = k_slope * loc(:, 2) + b_intercept;
t1wt2w_pred = t1wt2w + 0.5 * (y - t1wt2w);

a = 0.128 * t1wt2w_pred + 0.5549;
AMPA = NMDA./a - NMDA;
