subjects = {'sub-0002','sub-0003','sub-0004','sub-0006','sub-0007'};
base = 'D:\MEG-resting\MEG_raw\ds000247\ds000247';
data_idx = 28:297;
data_all_sub = cell(length(subjects),1);

idx = 1;  
subject = subjects{idx};

dataset_path = fullfile(base, subject, 'ses-01', 'meg', ...
    [subject '_ses-01_task-rest_run-01_meg.ds']);
hdr = ft_read_header(dataset_path);
data_labels = hdr.label(28:297);
data_meg= thedata; % 导入D:\MEG-resting文件夹下的S1.mat S2.mat S3.mat S4.mat S5.mat

fs = hdr.Fs;  

data_ft = [];
data_ft.trial{1} = data_meg;
data_ft.time{1}  = (0:size(data_meg,2)-1) / fs;
data_ft.label    = hdr.label(data_idx);
data_ft.fsample  = fs;

cfg_filt = [];
cfg_filt.bpfilter = 'yes';
cfg_filt.bpfreq   = [6 16];%[9 13];
cfg_filt.bpfiltord = 4;
data_filt = ft_preprocessing(cfg_filt, data_ft);
filt_data = data_filt.trial{1};  % 滤波之后的时间序列

window = 4 * fs;        % 4秒窗，频率分辨率 0.25 Hz
noverlap = window / 2;  % 50% 重叠
nfft = 2^nextpow2(window);

% 计算频率向量（取正频率部分）
freqs = (0:nfft/2-1) * (fs/nfft);

freq_range = freqs >= 6 & freqs <= 16;
freqs_plot = freqs(freq_range);

nCh = size(filt_data, 1);
pxx_all = zeros(nCh, sum(freq_range)); % PSD
for ch = 1:nCh
    [pxx, ~] = pwelch(filt_data(ch, :), window, noverlap, nfft, fs);
    pxx_all(ch, :) = pxx(freq_range)';
end