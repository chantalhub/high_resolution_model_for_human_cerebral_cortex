import os
import numpy as np
from scipy.sparse import load_npz, hstack, csr_matrix, save_npz, csc_matrix
import time

# 读取subject列表
subject_list_file = '/home/liuchang/High_resolution/structural/subject_list_HCP_255.txt'
with open(subject_list_file, 'r') as file:
    subject_list = [line.strip() for line in file]
row, col = load_npz('100206_retest_smooth_structural_connectivity.npz').shape
# 初始化一个存储矩阵的数组
matrix = load_npz('100206_retest_smooth_structural_connectivity.npz')
ratio = (matrix!=0).astype(int, copy=False)
# 循环处理每个subject
for subject in subject_list:
    print(subject)
    if subject == "100206":
        print("哈哈哈")
        pass
    else:
        # 构造文件路径
        npz_file = f'/home/liuchang/High_resolution/structural/{subject}_retest_smooth_structural_connectivity.npz'
        # 读取npz文件
        start = time.perf_counter()
        sc = load_npz(npz_file)
        end = time.perf_counter()
        print("load_npz: ", end-start)
        start = time.perf_counter()
        ratio = ratio + (sc!=0).astype(int, copy=False)
        end = time.perf_counter()
        print("ratio: ", end-start)
        start = time.perf_counter()
        matrix = matrix + sc
        end = time.perf_counter()
        print("matrix: ", end-start)
pos = (ratio/len(subject_list) > 0.1)
# 初始化组平均结构连接矩阵
group_average_sc = csr_matrix((row, col))
# 对于每个元素，如果大于50%的文件中该元素不等于0，取平均值；否则为0
group_average_sc[pos] = matrix[pos] / ratio[pos]
save_npz("group_average_sc_raw.npz", group_average_sc)
print("rows = {0}, cols = {1}".format(group_average_sc.shape[0], group_average_sc.shape[1]))
print("group_average_sc的密度：", group_average_sc.nnz / (group_average_sc.shape[0]*group_average_sc.shape[1]))



