v2：恢复“单邻居强度×邻居数量×平均方位权重”的社会输入标量

S_i=sum_j w_ij exp[cos(theta_ij-phi_i)-1]=n_i*bar(w_i)*C_i。
初始S0在遮挡集合冻结后、首次更新前计算，参考自身朝向取公共事件检查点。
A_i=sum_jw_ij作为总强度对照。完整高斯输入与神经动力学沿用父目录原模型运行器。
方案见社会输入三因素修订方案.v2.txt，参数见config.v2.json。

仅运行准备检查（模型推进步数为0）：
  python experiments/occlusion-input-consistency-20261030/v2/test_v2.py

未来获得明确运行指令后，在仓库根目录执行：
  addpath('experiments/occlusion-input-consistency-20261030');
  runOcclusionConsistency('model-smoke','model-smoke-v2');
  runOcclusionConsistency('execute','run-v2');
  python experiments/occlusion-input-consistency-20261030/v2/analyze_v2.py --run run-v2

原始运行结果及父目录v1文件保留；新的三因素表及图写入results/run-v2/analysis-v2/。
依赖沿用父目录requirements.txt。计算及绘图脚本与模型运行分离，分析脚本只读取原始轨迹与检查点。
主实验500个独立种子、两组历史对照各100个种子，共14000条计划轨迹；模型实验仍待运行。
当前分支YCwork，v2起点6fe205d。仅新增此目录，原模型、历史数据、论文与main保持原样。
本轮证据为定义修订和解析/构造数据软件检查。回退使用git revert本次提交。
