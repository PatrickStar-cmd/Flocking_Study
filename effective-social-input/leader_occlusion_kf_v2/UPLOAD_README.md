# 领航者方位缺测 KF 扩展实验 v2

这是 `10.13周/领航者实验` 的可审阅结果子集，上传到仓库时保留了代码、实验协议、汇总统计、主图和独立复核记录。

## 结果范围

- 仿真总数：858 条；独立逐运行复核：`results/checks/validation.json` 为 `PASS`。
- 主比较：N=20、3 个领航者、17 个跟随者、51 条领航者观测边缺测、遮挡 60 步。
- 评估种子：30 个；KF-mean 相对 Hold 的主场景 GO 差值为 `+0.0058506193`，95% 配对 t 区间为 `[+0.0054826421, +0.0062185965]`。
- 主结果图位于 `results/figures/`，逐种子配对数据位于 `results/analysis/`。

## 复现入口

从 `scripts/` 运行：

```powershell
python run_experiments.py --name repeat_v2 --workers 6
python verify_results.py results/repeat_v2
python analyze_results.py results/repeat_v2
```

具体参数和时间标签见 `docs/EXECUTION_PROTOCOL.md`、`docs/protocol_v2.md` 和 `docs/DATA_DICTIONARY.md`。

## 未上传的数据

逐步 HDF5/CSV、共享初态与完整测量噪声、183 页逐种子时间序列 PDF、183 页逐种子轨迹 PDF 以及大型时间序列 CSV 未上传。它们总计约 36 GB，保留在本地目录：

`C:\Users\WHQ\Studying\大三上学期\无人机群\10.13周\领航者实验\results\full_v2_run1`

这些文件不影响阅读汇总结果和复核状态；如需逐步数据复核，应从该本地目录取得对应文件。
