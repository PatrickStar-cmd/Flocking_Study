# P2.1 几何解耦与Q校准

正式报告：[P21_FINAL_REPORT.md](P21_FINAL_REPORT.md)。冻结停止判定的来源与边界：[DECISION_INTERPRETATION.md](DECISION_INTERPRETATION.md)。代码与运行方法：[实验代码README](../../../experiments/kf-p21-20261102/README.md)。

当前状态：A–E执行及审核完成，在全局Q选择门槛停止（STOP_NO_ADMISSIBLE_Q）。5个旧seed只作审计；20个合成噪声流、320个解析条件仅为机理/软件验证；12个新开发seed生成24条物理参考，保留1200个候选源窗口、7个不同配置共8400个方法/窗口条件。

独立核对504个开发原始预测文件SHA、全部58555个非空窗口RMSE、几何标签与纯真值三分位；覆盖判据另外从原始theta/P复算。q_star=null，原预留30个验证seed全部未使用。没有P3/P4或GO收益结论。

- baseline/：P2冻结来源SHA与270行复算。
- g0-r02/：正式位置标签与几何/NIS/尾部/覆盖分层；初版g0保留QA。正式全量尾部图在review/figures/G0_full_tails.png（原图样本截取仅作为绘图QA历史）。
- g1/：320个解析条件的完整真值、观测、缺失和三方法预测；正式机理年龄图见review/figures/G1_mechanism_age.png。
- development-r01/raw/input/：24条完整物理输入；labels/和GEOMETRY_FROZEN.json为开发纯真值分层。
- development-r01/raw/replay/：504个共享fresh及分支无损XOR文件。engine.recover_branch恢复所有候选的遮挡/恢复400步；预200步来自共享前缀。
- development-r01/analysis/：全部窗口、seed汇总、Q网格、年龄曲线、NIS尾部、P分箱和门限；review/为独立审核。
- development-r01/FROZEN_DECISION.json及.sha256：冻结参数、停止原因、未生成验证的声明及完整源码SHA。
- analysis/figures/：主要开发比较图与明确标注“未执行”的验证面板。
- logs/、QA_HISTORY.md、sources/：原始日志、QA版本、运行来源快照。

大型MAT/NPZ/压缩CSV采用Git LFS。原模型、P1/P2和历史结果只读，main未修改；本地结果准备完成不代表远端已发布，远端写入须遵照用户授权。公共回退使用对应提交git revert，保留历史。
