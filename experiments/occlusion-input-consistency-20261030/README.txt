遮挡后初始社会输入—最终运动一致性（v1，2026-10-30）

状态：方案与脚本准备，模型实验待运行。
方案：实验方案.v1.txt；参数：config.v1.json。

横轴：论文最终选用的q*=1，S0=mean_i sum_j w_ij=0.012*r，r=0..19。
纵轴：末1000步的mean_t min(GO(t),LO(t))；GO和LO均使用实际二维速度，LO包括自身和5个最近邻。
主实验500个新种子，20档遮挡，共10000条轨迹。额外两档历史对照各100种子，总计14000条。
全部模型过程使用原神经场动力学；缺失边删除；观测理想、领航线索为零。

MATLAB R2024b，在仓库根目录查看计划（零模型步）：
  addpath('experiments/occlusion-input-consistency-20261030');
  runOcclusionConsistency;

软件准备检查（静态分析及人工构造指标，零模型步）：
  checkPreparation;
  python experiments/occlusion-input-consistency-20261030/test_preparation.py

未来获得执行指令后，先进行模型回归smoke，再使用新的run_id：
  runOcclusionConsistency('model-smoke','model-smoke-v1');
  runOcclusionConsistency('execute','run-v1');
  python experiments/occlusion-input-consistency-20261030/analyze.py --run run-v1
以上execute命令属于未来实验入口。本次交付不执行该命令。
用seedOverride参数可单独生成一小批软件smoke轨迹；完整主分析要求全部冻结种子。

Python依赖：requirements.txt。每次运行复制源代码与参数到结果目录sources/。
结果目录results/、Python缓存及运行日志由本目录.gitignore隔离。
主模型脚本独立可移植，源公式与原始实现参考记录在manifest中。
本轮的静态检查及人工数据测试记录位于checks/，其中模型步数明确为0。

分支YCwork，起点6f1ee7b；仅上传此目录的方案与脚本。原始模型、历史结果、论文和main保持原样。
证据类型：预先设计的实验与软件准备检查。未来模型结果按实际完成状态另行登记。
回退：git revert对应上传提交，保留公共历史。
