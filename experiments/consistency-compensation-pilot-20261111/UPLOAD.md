# 仓库归档说明

本目录是 2026-11-11 完成的“小规模一致性实验-2”的独立归档，上传到 `WHQ_work` 分支。所有新增文件均位于本目录内。

## 包含内容

- MATLAB 仿真、方位预测、参数校准和一致性指标代码。
- Python 分析程序、固定依赖版本、实验配置及方案。
- 6 个开发种子在两个候选时序下的原始 MAT 轨迹、逐种子 CSV 和筛选记录。
- 24 个独立测试种子的原始 MAT 轨迹、逐种子 CSV、校准参数与完成清单。
- 配对统计、结果图、审计记录与结论说明。

上传时保留了原实验的文件结构与内容，仅排除 Python 缓存；另添加本说明、目录内 `.gitignore`、`.gitattributes` 和 `UPLOAD_INVENTORY.csv`。清单列出原实验 99 个归档文件的相对路径、字节数与 SHA256，可用于核对下载内容。目录内 `.gitattributes` 禁用换行符自动转换，确保下载后的文件字节与原实验一致。

## 查看和复核

先阅读 [README.md](README.md)，完整数字见 [测试集报告](results/test/analysis/README.md)。

在本目录安装 `requirements.txt` 中的 Python 依赖后，运行：

```powershell
python analyze_results.py test
```

该命令会读取本目录已有的全部测试轨迹，重新核对指标、遮挡时序和校准参数，并生成配对统计及图表；不需要重新运行 MATLAB 仿真。

重新仿真时请使用独立的新结果目录，参照 `PROTOCOL.md` 和 `README.md` 中的顺序执行，不要覆盖本次归档的校准参数和结果。

## 结论范围

本次是 20 架飞机的周期遮挡小规模实验。主条件下，数据支持 KF 相对 Drop 的一致性改善；尚不足以证明 KF 优于 Hold、Linear 或 Damped Holt，也没有验证通用的社会输入阈值。
