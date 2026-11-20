% 读取 subject_list_HCP.txt 文件
fid = fopen('/home/liuchang/High_resolution/structural/subject_list_HCP_255.txt', 'r');
subjects = textscan(fid, '%s');
fclose(fid);

subjects = subjects{1}; % 将 cell 转换为字符串数组

% 设置数据文件夹路径
dataFolder = '/home/liuchang/High_resolution/t1wt2w/data/';

% 初始化一个空矩阵来存储所有数据
t1wt2w = [];

% 迭代处理每个主体
for i = 1:length(subjects)
    subject = subjects{i};
    fprintf('Processing subject: %s\n', subject);

    % 构造文件路径
    filePath = fullfile(dataFolder, subject, [subject, '.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii']);
    
    % 读取 MyelinMap 数据
    nii = ft_read_cifti(filePath);
    data = nii.myelinmap_bc_msmall;

    % 将数据添加到 allData 矩阵中
    t1wt2w = [t1wt2w, data];
end

t1wt2w = median(t1wt2w, 2);
